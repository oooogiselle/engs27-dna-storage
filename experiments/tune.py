"""Compression tuning: variants of the context-mixing model, measured in bits/char."""
import glob, math, random, sys
import solution as S

texts = [open(f, "rb").read() for f in sorted(glob.glob("texts/*.txt"))]
corpus = b"\n\n".join(texts)
rng = random.Random(5)
excerpts = []
for n in [150, 250, 400, 700, 1000, 1500, 2500, 4000] * 2:
    a = rng.randrange(0, len(corpus) - n)
    excerpts.append(corpus[a:a + n])


class Var(S._Model):
    ORDERS = (1, 2, 3, 4, 6)
    LR0, LR1 = 0.02, 0.0       # lr = LR0 + LR1 / sqrt(1 + bytes)
    LIMIT = 255
    APMW = 0.5
    SEL = "node2"
    SPARSE = False
    WORD2 = False
    CONF = False

    def code_byte(self, coder, byte=None):
        h = self.hist
        cs = [("o0", 0)]
        for k in self.ORDERS:
            cs.append((k, bytes(h[-k:])) if len(h) >= k else (k, None))
        cs.append(("w", self.wordh))
        if self.WORD2:
            cs.append(("w2", (self.wordh, getattr(self, "prevword", 0))))
        if self.SPARSE:
            cs.append(("s", bytes(h[-3:-1])) if len(h) >= 3 else ("s", None))
        node = 1
        last = h[-1] if h else 0
        lr = self.LR0 + self.LR1 / math.sqrt(1 + len(h))
        for i in range(8):
            bit_in = None if byte is None else (byte >> (7 - i)) & 1
            slots, st = [], []
            for m, c in cs:
                if c is None:
                    slots.append(None); st.append(0.0)
                    if self.CONF: st.append(0.0)
                    continue
                key = (m, c, node)
                e = self.t.get(key)
                if e is None:
                    e = self.t[key] = [0.5, 0]
                slots.append(e)
                s = S._stretch(min(max(e[0], 1e-4), 1 - 1e-4)) if e[1] else 0.0
                st.append(s)
                if self.CONF:
                    st.append(s * e[1] / (e[1] + 3) if e[1] else 0.0)
            st.append(0.3)
            if self.SEL == "node2":
                sel = (node, min(len(h), 2))
            elif self.SEL == "node":
                sel = node
            else:  # node + how many orders have seen this context
                sel = (node, sum(1 for e in slots[1:] if e is not None and e[1] > 0))
            w = self.w.get(sel)
            if w is None:
                w = self.w[sel] = [0.3] * len(st)
            dot = sum(a * b for a, b in zip(w, st))
            p = S._squash(dot)
            bucket = 0 if dot <= -8 else (32 if dot >= 8 else int((dot + 8) * 2))
            ak = (last, node, bucket)
            a = self.apm.get(ak)
            pf = p if a is None else (1 - self.APMW) * p + self.APMW * a
            bit = coder.code(pf, bit_in)
            a0 = p if a is None else a
            self.apm[ak] = a0 + 0.05 * (bit - a0)
            err = bit - p
            for j in range(len(w)):
                w[j] += lr * err * st[j]
            for e in slots:
                if e is not None:
                    e[1] = min(e[1] + 1, self.LIMIT)
                    e[0] += (bit - e[0]) / (e[1] + 0.5)
            node = node * 2 + bit
        b = node & 0xFF
        h.append(b)
        if 65 <= b <= 90 or 97 <= b <= 122 or b >= 128:
            self.wordh = hash((self.wordh, b | 32))
        else:
            if self.wordh:
                self.prevword = self.wordh
            self.wordh = 0
        return b


def bits_for(data, cls):
    m, c = cls(), S._BitCoder()
    for b in data:
        c.code(m.p_more(), 1); m.more += 1
        m.code_byte(c, b)
    c.code(m.p_more(), 0)
    return len(c.finish())


def score(cls):
    # weight short/medium/long equally: mean bits/char over excerpts
    return sum(bits_for(t, cls) / len(t) for t in excerpts) / len(excerpts)


if __name__ == "__main__":
    variants = {
        "current": {},
        "lr .03+.3/sqrt": dict(LR0=0.01, LR1=0.3),
        "lr .005+.5/sqrt": dict(LR0=0.005, LR1=0.5),
        "sel orders-seen": dict(SEL="seen"),
        "sel node": dict(SEL="node"),
        "apm .25": dict(APMW=0.25),
        "conf inputs": dict(CONF=True),
        "sparse": dict(SPARSE=True),
        "word2": dict(WORD2=True),
        "orders 1-5": dict(ORDERS=(1, 2, 3, 4, 5)),
        "limit 30": dict(LIMIT=30),
    }
    for name, kw in variants.items():
        cls = type(name, (Var,), kw)
        print(f"{name:20s} {score(cls):.4f} bits/char", flush=True)
