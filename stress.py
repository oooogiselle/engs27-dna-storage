"""Stress test: many random excerpts (different strands each time) on a channel.
    python3 stress.py [real|sim] [N]"""
import random, sys
from bench import Channel, evaluate, texts
from sim import SimChannel

mode = sys.argv[1] if len(sys.argv) > 1 else "sim"
N = int(sys.argv[2]) if len(sys.argv) > 2 else 40
corpus = "\n\n".join(texts().values())
rng = random.Random(42)
res = []
for k in range(N):
    n = rng.choice([150, 250, 400, 700, 1000, 1500, 2500, 4000])
    a = rng.randrange(0, len(corpus) - n)
    t = corpus[a:a + n]
    ch = Channel() if mode == "real" else SimChannel(seed=1000 + k)
    res.append((n,) + evaluate(t, ch, name=f"excerpt{k}", verbose=False))
    s, d, e = res[-1][1:]
    if e > 0:
        print(f"  FAIL n={n} density {d:.2f} acc {1 - e:.3f}", flush=True)
by = {}
for n, s, d, e in res:
    by.setdefault(n, []).append(s)
print({n: round(sum(v) / len(v), 2) for n, v in sorted(by.items())})
print(f"{mode}: MEAN {sum(r[1] for r in res) / N:.3f}  failures {sum(r[3] > 0 for r in res)}/{N}")
