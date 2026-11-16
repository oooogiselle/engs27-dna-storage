"""Prototype: bitwise context-mixing compressor (lpaq-style, tiny).

Per bit, each model (order-k byte context, word context) predicts P(bit=1) from
adaptive counts; a logistic mixer combines them; an APM refines the output.
Encoder and decoder learn identically, so nothing but the bit stream is sent.
EOF is signalled by a 1-bit flag before each byte (cheap once learned).
"""
import math

def _stretch(p):
    return math.log(p / (1 - p))

def _squash(x):
    if x > 30: x = 30
    if x < -30: x = -30
    return 1 / (1 + math.exp(-x))

PREC = 1 << 16


class _BitCoder:
    """Binary arithmetic coder on 32-bit range, probabilities in 1/65536."""
    def __init__(self, bits=None):
        self.lo, self.hi = 0, 0xFFFFFFFF
        self.bits = bits
        if bits is None:
            self.out = []
        else:
            self.pos, self.x = 0, 0
            for _ in range(32):
                self.x = (self.x << 1) | self._nb()

    def _nb(self):
        b = self.bits[self.pos] if self.pos < len(self.bits) else 0
        self.pos += 1
        return b

    def code(self, p1, bit=None):
        """p1 = P(bit==1) as float. Encode bit, or decode and return it."""
        q = int(p1 * PREC)
        q = min(max(q, 1), PREC - 1)
        mid = self.lo + ((self.hi - self.lo) >> 16) * q + ((((self.hi - self.lo) & 0xFFFF) * q) >> 16)
        if bit is None:
            bit = 1 if self.x <= mid else 0
        if bit:
            self.hi = mid
        else:
            self.lo = mid + 1
        while (self.lo ^ self.hi) & 0x80000000 == 0:
            if self.bits is None:
                self.out.append(self.hi >> 31)
            else:
                self.x = ((self.x << 1) & 0xFFFFFFFF) | self._nb()
            self.lo = (self.lo << 1) & 0xFFFFFFFF
            self.hi = ((self.hi << 1) & 0xFFFFFFFF) | 1
        return bit

    def finish(self):
        # flush enough bits to pin a value in [lo, hi]
        for i in range(32):
            self.out.append((self.lo >> (31 - i)) & 1)
            # stop once prefix+zeros >= lo and prefix+ones <= hi
            pref = self.lo >> (31 - i) << (31 - i)
            if pref >= self.lo and (pref | ((1 << (31 - i)) - 1)) <= self.hi:
                break
        return self.out


class Model:
    def __init__(self, orders=(1, 2, 3, 4, 6), lr=0.02, limit=60, word=True):
        self.orders = orders
        self.word = word
        self.t = {}                     # (model, ctxhash, node) -> [p, n]
        self.n_in = len(orders) + 1 + (1 if word else 0) + 1   # + order0 + bias
        self.w = {}                     # mixer weights per selector
        self.lr, self.limit = lr, limit
        self.hist = []
        self.wordh = 0
        self.apm = {}

    def _ctxs(self):
        h = self.hist
        cs = [("o0", 0)]
        for k in self.orders:
            if len(h) >= k:
                cs.append((k, bytes(h[-k:])))
            else:
                cs.append((k, None))
        if self.word:
            cs.append(("w", self.wordh))
        return cs

    def code_byte(self, coder, byte=None):
        cs = self._ctxs()
        node = 1
        for i in range(8):
            bit_in = None if byte is None else (byte >> (7 - i)) & 1
            slots, st = [], []
            for m, c in cs:
                if c is None:
                    slots.append(None); st.append(0.0); continue
                key = (m, c, node)
                e = self.t.get(key)
                if e is None:
                    e = [0.5, 0]
                    self.t[key] = e
                slots.append(e)
                # confidence-scaled stretch: unseen nodes contribute 0
                st.append(_stretch(min(max(e[0], 1e-4), 1 - 1e-4)) if e[1] else 0.0)
            st.append(0.3)                                   # bias
            sel = (node, min(len(self.hist), 2))
            w = self.w.get(sel)
            if w is None:
                w = [0.3] * len(st)
                self.w[sel] = w
            dot = sum(a * b for a, b in zip(w, st))
            p = _squash(dot)
            # APM on order-1 context
            ak = (self.hist[-1] if self.hist else 0, node, int((dot + 8) * 2) if -8 < dot < 8 else (0 if dot <= -8 else 32))
            a = self.apm.get(ak)
            pf = p if a is None else 0.5 * p + 0.5 * a
            bit = coder.code(pf, bit_in)
            # update APM
            self.apm[ak] = (a if a is not None else p) + 0.05 * (bit - (a if a is not None else p))
            # update mixer
            err = bit - p
            for j in range(len(w)):
                w[j] += self.lr * err * st[j]
            # update models
            for e in slots:
                if e is None:
                    continue
                e[1] = min(e[1] + 1, self.limit)
                e[0] += (bit - e[0]) / (e[1] + 0.5)
            node = node * 2 + bit
        b = node & 0xFF
        self.hist.append(b)
        if (65 <= b <= 90) or (97 <= b <= 122) or b >= 128:
            self.wordh = hash((self.wordh, b | 32))
        else:
            self.wordh = 0
        return b


def compress(data, **kw):
    m, c = Model(**kw), _BitCoder()
    eof = Model(orders=(), word=False)      # tiny model for the "more bytes?" flag
    for b in data:
        c.code(_eof_p(eof), 1); _eof_up(eof, 1)
        m.code_byte(c, b)
    c.code(_eof_p(eof), 0)
    return c.finish()


def decompress(bits, limit=1 << 20, **kw):
    m, c = Model(**kw), _BitCoder(bits)
    eof = Model(orders=(), word=False)
    out = []
    while len(out) < limit:
        more = c.code(_eof_p(eof))
        _eof_up(eof, more)
        if not more:
            break
        out.append(m.code_byte(c))
    return bytes(out)


def _eof_p(m):
    n = m.t.get("n", 0)
    return (n + 1) / (n + 2)

def _eof_up(m, bit):
    m.t["n"] = m.t.get("n", 0) + bit


if __name__ == "__main__":
    import glob, itertools
    files = sorted(glob.glob("texts/*.txt"))
    texts = [open(f, "rb").read() for f in files]
    print("files", [f.split("/")[-1][:10] for f in files])
    for kw in [dict(), dict(lr=0.01), dict(lr=0.04), dict(limit=30), dict(limit=255),
               dict(orders=(1, 2, 3, 4)), dict(orders=(1, 2, 3, 5)), dict(word=False)]:
        row = []
        for t in texts:
            b = compress(t, **kw)
            row.append(f"{len(b) / len(t):.2f}")
        print(kw, " ".join(row), "avg", f"{sum(float(x) for x in row) / len(row):.3f}")
    t = texts[0]
    assert decompress(compress(t)) == t
    print("roundtrip ok")
