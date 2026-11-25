"""Ablation study: switch off one piece of the design at a time and measure what it was worth.

    python3 ablate.py sim [N]     # N excerpts per (variant, length) on the simulator (default 40)
    python3 ablate.py real [N]    # N excerpts per (variant, length) on the course server (default 5)

solution.py is not edited: each variant patches the module's globals inside a worker process
(its functions look globals up at call time, so encode/decode pick the patches up).
Every variant sees the same excerpts, so differences between rows are paired comparisons.
Writes results/ablation_<mode>.csv (one row per run) and results/ablation_<mode>_summary.txt.
"""
import contextlib, csv, importlib.util, math, os, random, sys, zlib
from multiprocessing import Pool

import solution as S
from bench import Channel, evaluate, texts
from sim import SimChannel

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CM = _load("solution_cm", "experiments/solution_cm.py")   # our earlier version: context mixing
FIXED = _load("partner_solution", "partner/solution.py")   # Huffman with a fixed English table
LENGTHS = [150, 250, 1000, 4000]
# Break rate per site for random (unconstrained) DNA, fitted on the real server from the fraction of
# full-length reads of 100 random strands at L = 24..60 (0.30, 0.18, 0.10, 0.086, 0.039): p ~ 0.053.
# (The simulator gives ~0.039 here, so it is optimistic about unconstrained strands.)
P_BREAK_PLAIN = 0.053


# ---------------------------------------------------------------- the planner's risk model
def fail_prob(L, n, r, nchars):
    """Planner's predicted P(passage unrecoverable): the XOR groups cannot rebuild every lost strand."""
    N = n + r
    lam = 4 * nchars / (N * (1 - S.P_LOST)) * (1 - S.P_BREAK) ** (L - 1)
    term, p_few = math.exp(-lam), 0.0
    for k in range(S.MIN_FULL):
        p_few += term
        term *= lam / (k + 1)
    q = S.P_LOST + (1 - S.P_LOST) * p_few
    return 1 - S._p_success(n, r, q)


def plan(nbits, nchars, r_max=S.MAX_R):
    """Same search as solution._plan, with the parity count capped at r_max."""
    best = None
    for L in S.L_RANGE:
        K = S._K(L)
        if K < 4:
            continue
        n = max(1, -(-nbits // (8 * K)))
        if n > 200:
            continue
        for r in [x for x in [0, 1] + list(range(3, S.MAX_R + 1)) if x <= r_max]:
            pf = fail_prob(L, n, r, nchars)
            if pf > S.MAX_FAIL:
                continue
            ev = (1 - pf) / ((n + r) * L)
            if best is None or ev > best[0]:
                best = (ev, L, n, r)
    return best[1:] if best else S._plan_fallback(nbits)


def fixed_plan(L, r):
    return lambda nbits, nchars: (L, max(1, -(-nbits // (8 * S._K(L)))), r)


# ---------------------------------------------------------------- replacement components
def _bytes_to_bits(data):
    return [bit for byte in data for bit in S._to_bits(byte, 8)]


def _bits_to_bytes(bits):
    return bytes(S._to_int(bits[i:i + 8]) for i in range(0, len(bits) // 8 * 8, 8))


def raw_compress(data):                      # no compression: UTF-8 bytes + a 0x00 terminator
    return _bytes_to_bits(data + b"\0")


def raw_decompress(bits, stop_bit=None, limit=1 << 20):
    data = _bits_to_bytes(bits if stop_bit is None else bits[:stop_bit])
    return data.split(b"\0", 1)[0][:limit]


def zlib_compress(data):                     # raw DEFLATE (no zlib header/checksum), level 9
    c = zlib.compressobj(9, zlib.DEFLATED, -15)
    return _bytes_to_bits(c.compress(data) + c.flush())


def zlib_decompress(bits, stop_bit=None, limit=1 << 20):
    data = _bits_to_bytes(bits if stop_bit is None else bits[:stop_bit])
    try:
        return zlib.decompressobj(-15).decompress(data, limit)
    except zlib.error:
        return b""


def fixed_compress(data):                    # Huffman from a built-in English table: no table sent
    return _bytes_to_bits(FIXED.compress(data.decode("utf-8")))


def fixed_decompress(bits, stop_bit=None, limit=1 << 20):
    return FIXED.decompress(_bits_to_bytes(bits if stop_bit is None else bits[:stop_bit])).encode("utf-8")[:limit]


def plain_capacity(L):                       # no constrained code: 2 bits per base
    return 2 * L


def plain_unrank(x, L):
    return "".join("ACGT"[(x >> 2 * (L - 1 - i)) & 3] for i in range(L))


def plain_rank(s):
    x = 0
    for ch in s:
        k = "ACGT".find(ch)
        if k < 0:
            return None
        x = 4 * x + k
    return x


VARIANTS = {
    "full": {},
    "raw-bytes": {"compress": raw_compress, "decompress": raw_decompress},
    "fixed-table": {"compress": fixed_compress, "decompress": fixed_decompress},
    "zlib": {"compress": zlib_compress, "decompress": zlib_decompress},
    "context-mixing": {"compress": CM.compress, "decompress": CM.decompress},
    "no-constraint": {"capacity": plain_capacity, "unrank": plain_unrank, "rank": plain_rank,
                      # short strands decode from 1-3 clean full reads, so the 8-read rule
                      # (tuned for L ~ 300) would be unfairly strict here
                      "P_BREAK": P_BREAK_PLAIN, "L_RANGE": range(20, 801, 4), "MIN_FULL": 3},
    "no-parity": {"_plan": lambda nb, nc: plan(nb, nc, r_max=0)},
    "fixed-plan": {"_plan": fixed_plan(200, 2)},
    "safe": {"MAX_FAIL": 0.01},
}
DESCRIBE = {
    "full": "everything on (the submission)",
    "raw-bytes": "no compression (UTF-8 bytes)",
    "fixed-table": "Huffman with a fixed English table (no table sent)",
    "zlib": "zlib/DEFLATE instead of Huffman",
    "context-mixing": "our earlier context-mixing coder instead of Huffman",
    "no-constraint": "plain 2 bits/base, GG/CC allowed",
    "no-parity": "no XOR parity strands",
    "fixed-plan": "fixed L=200, r=2 instead of per-passage plan",
    "safe": "risk cap MAX_FAIL=0.01",
}


@contextlib.contextmanager
def patched(variant):
    saved = {k: getattr(S, k) for k in VARIANTS[variant]}
    for k, v in VARIANTS[variant].items():
        setattr(S, k, v)
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(S, k, v)


# ---------------------------------------------------------------- running
_CORPUS = None


def excerpt(n, k):
    global _CORPUS
    if _CORPUS is None:
        _CORPUS = "\n\n".join(texts().values())
    a = random.Random(1000 * n + k).randrange(0, len(_CORPUS) - n)
    return _CORPUS[a:a + n]


def run(task):
    mode, variant, n, k = task
    text = excerpt(n, k)
    with patched(variant):
        bits = S.compress(text.encode("utf-8"))
        L, nd, r = S._plan(len(bits), len(text))
        pf = fail_prob(L, nd, r, len(text))
        ch = Channel() if mode == "real" else SimChannel(seed=100000 * n + k)
        score, density, e = evaluate(text, ch, verbose=False)
    return dict(variant=variant, chars=n, k=k, L=L, n_data=nd, r=r,
                bits_per_char=round(len(bits) / n, 3), p_fail_pred=round(pf, 4),
                density=round(density, 4), e=round(e, 4), score=round(score, 4))


def summarize(rows):
    lines = [f"{'variant':15s} " + " ".join(f"{n:>22d}" for n in LENGTHS) + "   mean",
             f"{'':15s} " + " ".join(f"{'score dens fail L/r':>22s}" for _ in LENGTHS)]
    for v in VARIANTS:
        cells, means = [], []
        for n in LENGTHS:
            rs = [x for x in rows if x["variant"] == v and x["chars"] == n]
            if not rs:
                cells.append(f"{'-':>22s}")
                continue
            sc = sum(x["score"] for x in rs) / len(rs)
            de = sum(x["density"] for x in rs) / len(rs)
            nf = sum(x["e"] > 0 for x in rs)
            Ls = sorted({x["L"] for x in rs})
            rr = sorted({x["r"] for x in rs})
            lr = f"{Ls[0]}" + (f"-{Ls[-1]}" if len(Ls) > 1 else "") + f"/{rr[0]}" + (f"-{rr[-1]}" if len(rr) > 1 else "")
            cells.append(f"{sc:5.2f} {de:4.2f} {nf:2d}/{len(rs):<3d}{lr:>7s}")
            means.append(sc)
        lines.append(f"{v:15s} " + " ".join(cells) + (f"  {sum(means) / len(means):5.2f}" if means else ""))
    lines.append("")
    lines.append("predicted vs observed failure rate (planner model):")
    for v in VARIANTS:
        rs = [x for x in rows if x["variant"] == v]
        if rs:
            pred = sum(x["p_fail_pred"] for x in rs) / len(rs)
            obs = sum(x["e"] > 0 for x in rs) / len(rs)
            lines.append(f"  {v:15s} predicted {pred:6.1%}   observed {obs:6.1%}  ({len(rs)} runs)")
    lines.append("")
    lines += [f"  {v:15s} {d}" for v, d in DESCRIBE.items()]
    return "\n".join(lines)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "sim"
    N = int(sys.argv[2]) if len(sys.argv) > 2 else (5 if mode == "real" else 40)
    # sanity: the copy of the planner matches solution._plan
    for nb, nc in [(500, 150), (4000, 1000), (15000, 4000)]:
        assert plan(nb, nc) == S._plan(nb, nc), (nb, nc)
    tasks = [(mode, v, n, k) for k in range(N) for n in LENGTHS for v in VARIANTS]
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    path = os.path.join(HERE, "results", f"ablation_{mode}.csv")
    rows = []
    with Pool(4 if mode == "real" else os.cpu_count()) as pool, open(path, "w", newline="") as f:
        w = None
        for i, row in enumerate(pool.imap_unordered(run, tasks), 1):
            if w is None:
                w = csv.DictWriter(f, fieldnames=list(row))
                w.writeheader()
            w.writerow(row)
            f.flush()
            rows.append(row)
            if i % 20 == 0 or i == len(tasks):
                print(f"{i}/{len(tasks)} runs", flush=True)
    rows.sort(key=lambda x: (list(VARIANTS).index(x["variant"]), x["chars"], x["k"]))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    text = summarize(rows)
    with open(os.path.join(HERE, "results", f"ablation_{mode}_summary.txt"), "w") as f:
        f.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
