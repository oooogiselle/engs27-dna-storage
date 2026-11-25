# ENGS 27 final project: storing text in DNA

Hand-in: **`solution.py`** (`encode(text) -> list[str]`, `decode(reads) -> str`) and the
one-page write-up in `writeup/`.

## Results (real course server)

| | density (bits/base) | accuracy | score |
|---|---|---|---|
| starter baseline (SAMPLE) | 1.44 | 0.92 | 0.87 |
| ours, 10 test passages (250–4000 chars) | 2.66 | 1.000 | **2.66** |
| ours, 20 random excerpts (150–4000 chars), ablation `full` row | 2.46 | 20/20 exact | 2.46 |
| ours, 24 random excerpts (150–4000 chars), notebook stress test | 2.50 | 23/24 exact | 2.39 |

Every technique is from the course (Huffman coding, parity checks and Hamming-style groups,
counting constrained strings, binomial/Poisson models). Failures, when they happen, are whole-passage
losses when more strands vanish than the XOR groups can rebuild (the channel drops ~3% of strands for
good): 5.6% of 160 simulated passages, mostly short ones that carry no parity on purpose. See "Risk knob" below.

## What each piece is worth (`python3 ablate.py real|sim`)

Real channel, mean score over 150 / 250 / 1000 / 4000-character excerpts (5 each), one piece changed at a time:

| variant | mean score | note |
|---|---|---|
| full design | **2.46** | 0/20 failed (5.6% in 160 simulated) |
| no compression | 1.50 | Huffman is worth 1.6× |
| fixed English Huffman table | 2.69 | no table sent: 2.91 vs 2.29 at 150 chars, same from 1000 on |
| GG/CC allowed, 2 bits/base | 0.20 | 17/20 failed: strands must be ~40 bases and many never arrive whole |
| no parity | 2.04 | same on short (planner uses none); 4/5 real, 28/40 simulated 4000-char passages fail |
| fixed L=200, r=2 | 1.41 | planner matters most for short passages (1.20 vs 2.29) |
| risk cap `MAX_FAIL=0.01` | 1.98 | simulated failures 5.6% → 0% |
| *context mixing (earlier version, not course material)* | 3.20 | what keeping to Huffman costs, mostly on long passages |

Full tables: `results/ablation_*_summary.txt`. Explanations for the check-in: `writeup/walkthrough.md`.

## What the channel does (Part 1 of the notebook)

* A strand breaks after **GG or CC** with p ≈ 0.25, elsewhere p ≈ 0.010 (independent).
* Reads shorter than 8 bases are dropped. ~0.6% errors per base (0.46% substitution, 0.11% insertion).
* ~3% of strands are **never read**. Which ones is fixed per sequence; no motif, end, GC or
  length effect was found.

## Design

compress (Huffman code built for the passage, codeword lengths sent first, ~4.5–6 bits/char)
→ fixed-length strands `[index | n data | r parity | payload | 8 parity bits]`
→ each strand unranked into a string with **no GG / CC** (1.83 bits/base)
→ r XOR parity strands, each covering a group of data strands; groups overlap like a Hamming code.
The encoder chooses strand length L and parity count r per passage to maximise
expected score. The decoder finds L from the read-length spike, clusters whole-strand reads,
takes a per-position majority, checks the parity bits, rebuilds lost strands group by group,
and Huffman-decodes.

## Files

| file | what |
|---|---|
| `solution.py` | the submission (self-contained, stdlib only) |
| `ENG_27_Project.ipynb` | finished notebook: channel exploration, solution cell, evaluation |
| `make_notebook.py` | builds the notebook from the starter + `solution.py` + `sim.py` |
| `sim.py` | channel simulator fitted to the measurements |
| `bench.py` | grader-style scoring: `python3 bench.py` (sim) or `python3 bench.py real` |
| `stress.py` | many random excerpts: `python3 stress.py real 40` |
| `ablate.py` | one-piece-off ablation study: `python3 ablate.py sim` (40/cell) or `python3 ablate.py real` (5/cell) |
| `partner/solution.py` | fixed-English-table Huffman prototype (used by the `fixed-table` ablation row) |
| `results/` | ablation CSVs (one row per run) and summary tables |
| `texts/` | test passages (`python3 texts/make_texts.py` regenerates them) |
| `experiments/` | prototypes (PPM vs context mixing, tuning), raw channel probes, and `solution_cm.py` (earlier context-mixing + Reed–Solomon version, used by the `context-mixing` ablation row) |
| `writeup/` | one-page write-up (HTML source + PDF) and `walkthrough.md` (check-in prep) |

## Workflow

Edit `solution.py`, then `python3 bench.py` and `python3 bench.py real`, then
`python3 make_notebook.py` to refresh the notebook's solution cell (this rebuilds the notebook
without outputs, including Part 4 from `results/ablation_real.csv`). Re-execute with
`jupyter nbconvert --to notebook --execute --inplace ENG_27_Project.ipynb` (needs numpy, matplotlib, requests).

## Risk knob

`MAX_FAIL` in `solution.py` caps the per-passage chance of an unrecoverable loss. At the
default (1.0) the encoder maximises expected score and short passages carry no parity
(~10% risk). Capping at 0.01 makes passages near-certain but lowers the expected score,
e.g. for 150-character passages from 2.29 to 1.65 on the real channel.

## Open question for the instructor

We assumed `solution.py` may not ship a pre-trained English model (everything the decoder
needs travels in the strands). If a fixed prior *is* allowed, short passages would compress
much better. Worth asking.
