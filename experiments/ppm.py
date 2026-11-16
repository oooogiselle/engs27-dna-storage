"""Prototype: adaptive PPM (method D, with exclusion) + binary arithmetic coder.

Nothing is transmitted except the bit stream: encoder and decoder build the same
model as they go. Symbols are bytes 0..255 plus EOF (256), so no length header.
"""

EOF = 256
NSYM = 257

# order -1 prior: printable ASCII + newline/tab are likely, other bytes rare.
# Generic "this is text" knowledge only, no English statistics.
_PRIOR = [0] * NSYM
for _s in range(NSYM):
    if 32 <= _s < 127 or _s in (9, 10):
        _PRIOR[_s] = 16
    elif _s == EOF:
        _PRIOR[_s] = 1
    else:
        _PRIOR[_s] = 1


class _Coder:
    """32-bit arithmetic coder (Witten-Neal-Cleary), bit-level output."""
    TOP = (1 << 32) - 1
    HALF = 1 << 31
    QTR = 1 << 30

    def __init__(self, bits=None):
        self.lo, self.hi = 0, self.TOP
        if bits is None:                       # encoder
            self.out, self.pending = [], 0
        else:                                  # decoder
            self.bits, self.pos = bits, 0
            self.val = 0
            for _ in range(32):
                self.val = (self.val << 1) | self._next()

    def _next(self):
        b = self.bits[self.pos] if self.pos < len(self.bits) else 0
        self.pos += 1
        return b

    def _emit(self, b):
        self.out.append(b)
        self.out.extend([1 - b] * self.pending)
        self.pending = 0

    def encode(self, cl, ch, tot):
        r = self.hi - self.lo + 1
        self.hi = self.lo + r * ch // tot - 1
        self.lo = self.lo + r * cl // tot
        while True:
            if self.hi < self.HALF:
                self._emit(0)
            elif self.lo >= self.HALF:
                self._emit(1); self.lo -= self.HALF; self.hi -= self.HALF
            elif self.lo >= self.QTR and self.hi < 3 * self.QTR:
                self.pending += 1; self.lo -= self.QTR; self.hi -= self.QTR
            else:
                break
            self.lo <<= 1; self.hi = (self.hi << 1) | 1

    def finish(self):
        self.pending += 1
        self._emit(0 if self.lo < self.QTR else 1)
        return self.out

    def target(self, tot):
        r = self.hi - self.lo + 1
        return ((self.val - self.lo + 1) * tot - 1) // r

    def consume(self, cl, ch, tot):
        r = self.hi - self.lo + 1
        self.hi = self.lo + r * ch // tot - 1
        self.lo = self.lo + r * cl // tot
        while True:
            if self.hi < self.HALF:
                pass
            elif self.lo >= self.HALF:
                self.lo -= self.HALF; self.hi -= self.HALF; self.val -= self.HALF
            elif self.lo >= self.QTR and self.hi < 3 * self.QTR:
                self.lo -= self.QTR; self.hi -= self.QTR; self.val -= self.QTR
            else:
                break
            self.lo <<= 1; self.hi = (self.hi << 1) | 1
            self.val = (self.val << 1) | self._next()


class _PPM:
    def __init__(self, order):
        self.order = order
        self.ctx = {}                          # context bytes -> {sym: count}

    def _dists(self, hist):
        """Yield (context-dict) from highest order to 0."""
        for k in range(min(self.order, len(hist)), -1, -1):
            yield k, self.ctx.get(bytes(hist[len(hist) - k:]) if k else b"")

    def _update(self, hist, s):
        for k in range(min(self.order, len(hist)), -1, -1):
            key = bytes(hist[len(hist) - k:]) if k else b""
            d = self.ctx.setdefault(key, {})
            d[s] = d.get(s, 0) + 1

    # Each step: a list of (symbol -> freq) tables to try in order, with escape.
    def code(self, coder, hist, s=None):
        excl = set()
        for k, d in self._dists(hist):
            if not d:
                continue
            items = [(x, 2 * c - 1) for x, c in d.items() if x not in excl]
            if not items:
                continue
            esc = len(items)
            tot = sum(f for _, f in items) + esc
            if s is not None:                  # encoding
                cl = 0
                for x, f in items:
                    if x == s:
                        coder.encode(cl, cl + f, tot)
                        return s
                    cl += f
                coder.encode(cl, cl + esc, tot)
            else:                              # decoding
                t = coder.target(tot)
                cl = 0
                for x, f in items:
                    if t < cl + f:
                        coder.consume(cl, cl + f, tot)
                        return x
                    cl += f
                coder.consume(cl, cl + esc, tot)
            excl.update(x for x, _ in items)
        # order -1
        items = [(x, _PRIOR[x]) for x in range(NSYM) if x not in excl]
        tot = sum(f for _, f in items)
        cl = 0
        if s is not None:
            for x, f in items:
                if x == s:
                    coder.encode(cl, cl + f, tot)
                    return s
                cl += f
        t = coder.target(tot)
        for x, f in items:
            if t < cl + f:
                coder.consume(cl, cl + f, tot)
                return x
            cl += f
        raise ValueError("bad stream")


def compress(data, order=4):
    m, c, hist = _PPM(order), _Coder(), []
    for s in list(data) + [EOF]:
        m.code(c, hist, s)
        m._update(hist, s)
        hist.append(s)
    return c.finish()


def decompress(bits, order=4, limit=1 << 20):
    m, c, hist, out = _PPM(order), _Coder(bits), [], []
    while len(out) < limit:
        s = m.code(c, hist)
        if s == EOF:
            break
        m._update(hist, s)
        hist.append(s)
        out.append(s)
    return bytes(out)


if __name__ == "__main__":
    import glob, sys
    for order in (2, 3, 4, 5):
        row = []
        for f in sorted(glob.glob("texts/*.txt")):
            t = open(f, "rb").read()
            b = compress(t, order)
            assert decompress(b, order) == t, f
            row.append(f"{len(b) / len(t):.2f}")
        print("order", order, " ".join(row))
    print("files", [f.split("/")[-1][:10] for f in sorted(glob.glob("texts/*.txt"))])
