"""Build ENG_27_Project.ipynb from the starter notebook + solution.py + sim.py.

The solution cell is generated from solution.py so the two can never drift apart.
    python3 make_notebook.py
then execute it:
    jupyter nbconvert --to notebook --execute --inplace ENG_27_Project.ipynb
"""
import copy, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
nb = json.load(open(os.path.join(HERE, "ENG_27_Project_Starter.ipynb")))
starter = nb["cells"]


def md(s):
    return {"cell_type": "markdown", "metadata": {}, "source": s.strip("\n")}


def code(s):
    return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [], "source": s.strip("\n")}


def src(i):
    c = starter[i]
    return "".join(c["source"]) if isinstance(c["source"], list) else c["source"]


cells = []
cells += [copy.deepcopy(c) for c in starter[0:2]]
# setup: keep the starter cell, add a figures folder and a fixed seed
cells.append(code(src(2) + "\nimport os\nos.makedirs('figures', exist_ok=True)\n"))
cells += [copy.deepcopy(c) for c in starter[3:5]]
for c in cells:
    c["outputs"] = [] if c["cell_type"] == "code" else c.get("outputs")
    if c["cell_type"] == "code":
        c["execution_count"] = None
    else:
        c.pop("outputs", None)

cells.append(md("""
## Part 1: Explore the channel

Everything below talks to the real server. Each experiment sends strands we know and
compares the reads with them. Results are summarised at the end of Part 1.
"""))
cells.append(code(src(6)))
cells.append(md(src(7)))
cells.append(code(src(8)))

cells.append(md("""
### 1.1 Are breaks independent? Read lengths vs. a geometric fit

If every site broke independently with one probability $p$, piece lengths would be geometric
(a straight line on a log scale). Reads shorter than 8 bases never come back.
"""))
cells.append(code("""
counts = collections.Counter(lengths)
xs = np.arange(8, 121)
ys = np.array([counts.get(x, 0) for x in xs])
p_hat = 1 / (np.mean(lengths) - 7)                       # geometric fit above the 8-base floor
fit = len(lengths) * p_hat * (1 - p_hat) ** (xs - 8)
plt.semilogy(xs, np.maximum(ys, 0.5), ".", label="real channel")
plt.semilogy(xs, fit, "-", label=f"geometric, p = {p_hat:.3f}")
plt.xlabel("read length"); plt.ylabel("count (log)"); plt.legend(); plt.title("Random 120-base strands")
plt.savefig("figures/read_lengths.png", dpi=150, bbox_inches="tight"); plt.show()
print(f"shortest read {min(lengths)}, reads of full length 120: {counts.get(120, 0)}")
"""))

cells.append(md("""
### 1.2 Does the chance of a break depend on the bases around it?

For every located read, the end of the read is a cut. Compare how often each local context
appears at a cut with how often it appears in the strands. A ratio of 1 means "no effect".
"""))
cells.append(code("""
rng = random.Random(3)
ctx_strands = ["".join(rng.choice("ACGT") for _ in range(200)) for _ in range(100)]
ctx_reads = channel.reads(ctx_strands, 40000)
ctx_loc = locate(ctx_reads, ctx_strands)
print(f"located {len(ctx_loc)} of {len(ctx_reads)} reads")

def cut_ratio(left, right, strands=ctx_strands, loc=ctx_loc, min_occ=30):
    L = len(strands[0])
    cuts, occ = collections.Counter(), collections.Counter()
    for j, p, n in loc:
        e = p + n
        if left <= e <= L - right and e < L:
            cuts[strands[j][e - left:e + right]] += 1
    for s in strands:
        for e in range(left, L - right):
            occ[s[e - left:e + right]] += 1
    tc, to = sum(cuts.values()), sum(occ.values())
    return {k: (cuts[k] / tc) / (occ[k] / to) for k in occ if occ[k] >= min_occ}

for left, right in [(1, 0), (0, 1), (2, 0), (0, 2), (3, 0), (2, 2)]:
    r = cut_ratio(left, right)
    srt = sorted(r.items(), key=lambda kv: -kv[1])
    print(f"{left} left/{right} right: ratio range {min(r.values()):.2f}-{max(r.values()):.2f}  "
          f"top {[(k, round(v, 1)) for k, v in srt[:4]]}")
"""))
cells.append(code("""
r2 = cut_ratio(2, 0)
keys = sorted(r2, key=lambda k: -r2[k])
plt.figure(figsize=(8, 3))
plt.bar(keys, [r2[k] for k in keys], color=["C3" if k in ("GG", "CC") else "C0" for k in keys])
plt.axhline(1, color="k", lw=0.5)
plt.ylabel("cut rate / average"); plt.xlabel("two bases just before the cut")
plt.title("Breaks happen after GG and CC")
plt.savefig("figures/break_context.png", dpi=150, bbox_inches="tight"); plt.show()
"""))
cells.append(md("""
The two bases **left** of the cut decide almost everything: after `GG` or `CC` the cut rate is
about 6x average, everywhere else it is about 0.25x. The bases to the right of the cut do not
matter (ratios stay near 1), and a third base on the left adds nothing beyond the dinucleotide.
So a long random strand is effectively chopped at GG/CC sites, which is why reads average
~25 bases. A strand that never contains `GG` or `CC` should break far less.
"""))

cells.append(md("""
### 1.3 How often does a read cover a whole strand of length $L$?

Compare random strands with strands that never contain `GG`/`CC`. If the remaining breaks are
independent with probability $p$ per site, a read is the whole strand with probability
$(1-p)^{L-1}$.
"""))
cells.append(code("""
def constrained_strand(L, rng):
    s = [rng.choice("ACGT")]
    while len(s) < L:
        b = rng.choice("ACGT")
        if not (b == s[-1] and b in "GC"):
            s.append(b)
    return "".join(s)

rng = random.Random(5)
Ls = [30, 60, 120, 200, 300]
full_rand, full_ban = [], []
for L in Ls:
    for ban, out in ((False, full_rand), (True, full_ban)):
        S = [constrained_strand(L, rng) if ban else "".join(rng.choice("ACGT") for _ in range(L)) for _ in range(60)]
        R = channel.reads(S, 10000)
        out.append(sum(len(r) == L for r in R) / len(R))
p_fit = np.mean([1 - f ** (1 / (L - 1)) for L, f in zip(Ls, full_ban)])
print("full-length fraction, random     :", [round(f, 3) for f in full_rand])
print("full-length fraction, no GG/CC   :", [round(f, 3) for f in full_ban])
print(f"fitted p per site with no GG/CC = {p_fit:.4f}")
xs = np.arange(20, 320)
plt.plot(Ls, full_rand, "o-", label="random strands")
plt.plot(Ls, full_ban, "s-", label="no GG / CC")
plt.plot(xs, (1 - p_fit) ** (xs - 1), "k--", lw=1, label=f"$(1-p)^{{L-1}}$, p={p_fit:.4f}")
plt.xlabel("strand length L"); plt.ylabel("fraction of reads = whole strand"); plt.legend()
plt.savefig("figures/full_coverage.png", dpi=150, bbox_inches="tight"); plt.show()
"""))

cells.append(md("""
### 1.4 Errors inside reads

About one read in ten is not an exact substring of any strand. Align those reads to their
source to see what kind of errors they contain.
"""))
cells.append(code("""
def semi_global(r, s):
    \"\"\"Align all of read r to any part of s. Returns counts of subs / insertions / deletions.\"\"\"
    n, m = len(r), len(s)
    D = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        D[i][0] = i
        for j in range(1, m + 1):
            D[i][j] = min(D[i-1][j] + 1, D[i][j-1] + 1, D[i-1][j-1] + (r[i-1] != s[j-1]))
    j = min(range(m + 1), key=lambda j: D[n][j]); i = n
    ops = collections.Counter()
    while i > 0:
        if j > 0 and D[i][j] == D[i-1][j-1] + (r[i-1] != s[j-1]):
            ops["sub" if r[i-1] != s[j-1] else "ok"] += 1; i -= 1; j -= 1
        elif D[i][j] == D[i-1][j] + 1:
            ops["ins"] += 1; i -= 1
        else:
            ops["del"] += 1; j -= 1
    return ops

K = 7
seed_index = collections.defaultdict(list)
for j, s in enumerate(ctx_strands):
    for i in range(len(s) - K + 1):
        seed_index[s[i:i+K]].append((j, i))
sample = ctx_reads[:3000]
is_exact = [any(r in s for s in ctx_strands) for r in sample]
bad = [r for r, ok in zip(sample, is_exact) if not ok]
ops = collections.Counter()
for r in bad:
    votes = collections.Counter((j, p - i) for i in range(len(r) - K + 1) for j, p in seed_index.get(r[i:i+K], []))
    if votes:
        (j, off), _ = votes.most_common(1)[0]
        ops += semi_global(r, ctx_strands[j][max(0, off - 4):off + len(r) + 6])
nall = sum(len(r) for r in sample)
print(f"reads with an error: {1 - sum(is_exact) / len(sample):.3f}")
print(f"per base over all reads: substitution {ops['sub'] / nall:.4f}, insertion {ops['ins'] / nall:.4f}, deletion {ops['del'] / nall:.5f}")

by_len = collections.defaultdict(lambda: [0, 0])
for r, ok in zip(sample, is_exact):
    b = min(len(r) // 20, 5); by_len[b][0] += 1; by_len[b][1] += not ok
print("error rate by read length:", {f"{20*b}-{20*b+19}": round(v[1] / v[0], 3) for b, v in sorted(by_len.items())})
"""))

cells.append(md("""
### 1.5 Some strands are never read at all

Send 900 strands and ask for 20,000 reads (~22 per strand). Count strands that get no reads,
resend the lost ones to check whether it is the *sequence* that is lost, and repeat at
several lengths.
"""))
cells.append(code("""
def reads_per_strand(S, R, k=12):
    idx = collections.defaultdict(set)
    for j, s in enumerate(S):
        for i in range(len(s) - k + 1):
            idx[s[i:i+k]].add(j)
    c = collections.Counter()
    for r in R:
        h = idx.get(r[:k], set()) if len(r) >= k else set()
        if len(h) == 1:
            c[next(iter(h))] += 1
    return [c[j] for j in range(len(S))]

rng = random.Random(11)
S = [constrained_strand(100, rng) for _ in range(900)]
cnt = reads_per_strand(S, channel.reads(S, 20000))
lost = [S[j] for j in range(len(S)) if cnt[j] == 0]
kept = [S[j] for j in range(len(S)) if cnt[j] >= 15][:200]
print(f"strands with zero reads: {len(lost)} of {len(S)} ({len(lost) / len(S):.1%})")
again = reads_per_strand(lost + kept, channel.reads(lost + kept, 10000))
print(f"resent: previously lost still lost {sum(c == 0 for c in again[:len(lost)])}/{len(lost)}, "
      f"previously kept now lost {sum(c == 0 for c in again[len(lost):])}/{len(kept)}")
loss_by_L = {}
for L in (40, 100, 200, 400):
    n = min(800, 90000 // L)
    S = [constrained_strand(L, rng) for _ in range(n)]
    c = reads_per_strand(S, channel.reads(S, 20 * n))
    loss_by_L[L] = sum(x == 0 for x in c) / n
print("loss rate by strand length:", {L: f"{v:.1%}" for L, v in loss_by_L.items()})
"""))
cells.append(md("""
### Summary of the channel

| what | measured | consequence for the design |
|---|---|---|
| breaks | after `GG` or `CC`: p ≈ 0.25; elsewhere p ≈ 0.010, independent | never write `GG`/`CC` (costs 8.6%: 1.83 instead of 2 bits/base) |
| full reads | with no GG/CC: $(1-p)^{L-1}$ with p ≈ 0.010 (74% at L=30, 30% at L=120) | strands of 250–450 bases still get dozens of whole-strand reads |
| read errors | ≈0.5% per base: ~85% substitutions, ~15% insertions | vote across many whole-strand reads; CRC on every strand |
| lost strands | ≈3% of strands never come back, fixed per sequence, not tied to length or motif | Reed–Solomon parity strands across strands |
| coverage | 4 reads/char at ~2 bases/char → ~2 reads per stored base, each 50–100 bases long | coverage is plentiful: compression is the main lever |
"""))

cells.append(md("""
## Part 2: Our solution

* **Compression.** A small bitwise context-mixing model (orders 0–4 and 6, plus a word
  model, mixed by an online logistic mixer, with an APM stage), driving a binary arithmetic coder.
  It is adaptive: the decoder learns the same model while decoding, so no table is sent.
  An end-of-text flag replaces a length header. About 3.5–4 bits/char on English, against 8 for raw bytes.
* **Constrained code.** Each strand's bits are one big integer, *unranked* into the
  set of length-$L$ strings with no `GG` and no `CC` (enumerative coding, exactly
  $\\lfloor\\log_2 N(L)\\rfloor \\approx 1.83L$ bits). One strand error never spreads to another strand.
* **Strand format.** `[index 8 bits | number of data strands 8 bits | payload | CRC-8]`.
* **Erasure code.** Systematic Cauchy Reed–Solomon over GF(256) across strands: any `r` lost
  strands can be rebuilt from `r` parity strands.
* **Choosing L and r.** For each passage the encoder picks the strand length and number of parity
  strands that maximise expected score, $\\text{density} \\times P(\\text{at most } r \\text{ strands fail})$,
  where a strand fails if it is lost (3.5%) or gets fewer than 8 whole-strand reads (Poisson model).
* **Decoding.** The strand length shows up as a spike in the read-length histogram. Whole-strand
  reads are clustered by Hamming distance, each cluster gets a per-position majority vote, the
  consensus is ranked back to bits and CRC-checked, missing strands are rebuilt from parity, and
  the stream is decompressed. If too many strands are missing, we decode up to the first gap.

**This cell is the solution.** Running it saves `solution.py`, the file we hand in.
"""))
cells.append(code("%%writefile solution.py\n" + open(os.path.join(HERE, "solution.py")).read()))

cells.append(md("""
## Part 3: Test it the way it will be graded

`LocalChannel` breaks every site with p = 0.03 and so is much harsher on long strands than the
real channel. We test on a simulator fitted to the Part 1 measurements, then on the real server,
across passages of different lengths and styles.
"""))
cells.append(code(open(os.path.join(HERE, "sim.py")).read().split('"""', 2)[2].strip()))
cells.append(code(src(12).replace("evaluate(SAMPLE, LocalChannel(seed=1))\nevaluate(SAMPLE, Channel())", "")
                  .rstrip() + "\n\nevaluate(SAMPLE, SimChannel(seed=1))\nevaluate(SAMPLE, Channel())"))
cells.append(code("""
import glob
passages = {os.path.basename(f)[:-4]: open(f, encoding="utf-8").read() for f in sorted(glob.glob("texts/*.txt"))}
results = {}
for name, text in passages.items():
    print(f"{name:18s} {len(text):5d} chars: ", end="")
    results[name] = evaluate(text, Channel())
print(f"\\nmean score on the real channel: {np.mean(list(results.values())):.3f}  (baseline: 0.87 on SAMPLE)")
"""))
cells.append(code("""
# Stress test: random excerpts mean different strands, so different strands get lost.
rng = random.Random(42)
corpus = "\\n\\n".join(passages.values())
stress = collections.defaultdict(list)
for k in range(24):
    n = rng.choice([150, 400, 1000, 2500, 4000])
    a = rng.randrange(0, len(corpus) - n)
    print(f"excerpt {k:2d} ({n} chars): ", end="")
    stress[n].append(evaluate(corpus[a:a + n], Channel()))
print({n: round(np.mean(v), 2) for n, v in sorted(stress.items())},
      f"failures: {sum(s == 0 for v in stress.values() for s in v)}/24")
"""))
cells.append(md(src(13)))

nb["cells"] = cells
json.dump(nb, open(os.path.join(HERE, "ENG_27_Project.ipynb"), "w"), indent=1)
print("wrote ENG_27_Project.ipynb with", len(cells), "cells")
