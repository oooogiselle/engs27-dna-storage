"""Score solution.py the way the grader does, on the simulator or the real server.

    python3 bench.py              # simulator, all texts, 3 seeds
    python3 bench.py real         # real course server, all texts, 1 run each
"""
import glob, importlib, math, os, sys, time

import solution
from sim import SimChannel

HERE = os.path.dirname(os.path.abspath(__file__))
SERVER_URL = "https://dna.cedar.host.dartmouth.edu"


class Channel:
    def __init__(self, url=SERVER_URL, key=""):
        import requests
        self.requests = requests
        self.url, self.key = url, key
        self.max_per_request = requests.get(url + "/api/health", timeout=30).json()["limits"]["max_reads_per_request"]

    def reads(self, strands, n):
        got = []
        while len(got) < n:
            ask = min(n - len(got), self.max_per_request)
            r = self.requests.post(self.url + "/api/reads", timeout=300,
                                   json={"key": self.key, "strands": list(strands), "n": ask})
            if r.status_code != 200:
                raise RuntimeError(r.text)
            got += r.json()["reads"]
        return got


def edit_distance(a, b):
    if a == b:
        return 0
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def H_info(e):
    return 0.0 if e >= 0.5 else (1.0 if e == 0 else 1 + e * math.log2(e) + (1 - e) * math.log2(1 - e))


def evaluate(text, channel, reads_per_char=4, verbose=True, name=""):
    strands = solution.encode(text)
    bases = sum(map(len, strands))
    reads = channel.reads(strands, reads_per_char * len(text))
    decoded = solution.decode(reads)
    density = 8 * len(text.encode()) / bases
    e = min(1.0, edit_distance(decoded[:2 * len(text)], text) / len(text))
    score = density * H_info(e)
    if verbose:
        print(f"{name:18s} {len(text):5d} ch {len(strands):3d} x {len(strands[0]):3d} = {bases:6d} bases  "
              f"density {density:.3f}  acc {1 - e:.4f}  SCORE {score:.3f}", flush=True)
    return score, density, e


def texts():
    return {os.path.basename(f)[:-4]: open(f, encoding="utf-8").read()
            for f in sorted(glob.glob(os.path.join(HERE, "texts", "*.txt")))}


if __name__ == "__main__":
    real = len(sys.argv) > 1 and sys.argv[1] == "real"
    seeds = [0] if real else [0, 1, 2]
    res = []
    t0 = time.time()
    for seed in seeds:
        for name, t in texts().items():
            ch = Channel() if real else SimChannel(seed=seed)
            res.append(evaluate(t, ch, name=name))
    n = len(res)
    print(f"MEAN score {sum(r[0] for r in res) / n:.3f}  density {sum(r[1] for r in res) / n:.3f}  "
          f"failures(e>0) {sum(r[2] > 0 for r in res)}/{n}  ({time.time() - t0:.0f}s)")
