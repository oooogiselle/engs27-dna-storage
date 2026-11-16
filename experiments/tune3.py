import solution as S, math
from tune import Var, excerpts
class M(Var):
    LR0, LR1 = 0.01, 0.3
    MINLEN = 5
    def __init__(self):
        super().__init__(); self.mtab = {}; self.mptr = None; self.mlen = 0
    def code_byte(self, coder, byte=None):
        h = self.hist
        # predicted byte from match model
        pred = h[self.mptr] if self.mptr is not None and self.mptr < len(h) else None
        self._pred = pred
        return Var.code_byte(self, coder, byte)
# patch: inject match input by overriding stretch list creation is easiest by copying code; do it via a wrapper on _squash inputs
import tune
src = open(__import__("os").path.join(__import__("os").path.dirname(__file__), "tune.py")).read().split("if __name__")[0]
src = src.replace("            st.append(0.3)\n", """            st.append(0.3)
            pr = getattr(self, '_pred', None)
            if pr is not None and ((pr | 256) >> (8 - i)) == node:
                bitp = (pr >> (7 - i)) & 1
                s = min(self.mlen, 24) / 4.0
                st.append(s if bitp else -s)
            else:
                st.append(0.0)
""").replace("            sel = (node, min(len(h), 2))", "            sel = (node, min(len(h), 2), pr is not None and ((pr | 256) >> (8 - i)) == node and min(self.mlen, 16) // 4)")
ns = {}; exec(compile(src, "t", "exec"), ns)
Base = ns["Var"]
class MM(Base):
    LR0, LR1 = 0.01, 0.3
    MINLEN = 5
    def __init__(self):
        super().__init__(); self.mtab = {}; self.mptr = None; self.mlen = 0
    def code_byte(self, coder, byte=None):
        h = self.hist
        self._pred = h[self.mptr] if self.mptr is not None and self.mptr < len(h) else None
        b = Base.code_byte(self, coder, byte)
        # update match state
        if self.mptr is not None and self._pred == b:
            self.mlen += 1; self.mptr += 1
        else:
            self.mlen = 0; self.mptr = None
        k = self.MINLEN
        if len(h) >= k:
            key = bytes(h[-k:])
            if self.mptr is None and key in self.mtab:
                self.mptr = self.mtab[key]; self.mlen = 0
            self.mtab[key] = len(h)
        return b
for ml in (4, 5, 6, 8):
    cls = type("m", (MM,), dict(MINLEN=ml))
    print("match minlen", ml, round(ns["score"](cls), 4), flush=True)
