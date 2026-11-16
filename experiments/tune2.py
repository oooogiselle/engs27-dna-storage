import tune, solution as S
from tune import Var, score
class V2(Var):
    LR0, LR1 = 0.01, 0.3
    W0 = 0.3; D = 0.5; ARATE = 0.05
    def code_byte(self, coder, byte=None):
        # reuse Var, but patch initial weights and counter rate via attributes
        return Var.code_byte(self, coder, byte)
# monkeypatch helpers inside Var through closures is messy; vary via subclass attributes:
import types
src = open(__import__("os").path.join(__import__("os").path.dirname(__file__), "tune.py")).read()
src = src.replace("[0.3] * len(st)", "[self.W0] * len(st)").replace("(e[1] + 0.5)", "(e[1] + self.D)").replace("0.05 * (bit - a0)", "self.ARATE * (bit - a0)")
src = src.split('if __name__')[0].replace("class Var(S._Model):", "class Var(S._Model):\n    W0 = 0.3; D = 0.5; ARATE = 0.05")
ns = {}; exec(compile(src, "tunep", "exec"), ns)
Base = ns["Var"]; sc = ns["score"]
for name, kw in {"base(lrdecay)": {}, "w0 .15": dict(W0=0.15), "w0 .2": dict(W0=0.2), "w0 .45": dict(W0=0.45),
                 "D 1.0": dict(D=1.0), "D 0.3": dict(D=0.3), "D 1.5": dict(D=1.5),
                 "arate .02": dict(ARATE=0.02), "arate .1": dict(ARATE=0.1), "lr .01+.5": dict(LR1=0.5)}.items():
    kw = dict(LR0=0.01, LR1=0.3, **kw) if "LR1" not in kw else dict(LR0=0.01, **kw)
    print(f"{name:15s} {sc(type(name,(Base,),kw)):.4f}", flush=True)
