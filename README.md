# ENGS 27 final project: storing text in DNA

Hand-in: **`solution.py`** (`encode(text) -> list[str]`, `decode(reads) -> str`) and the
one-page write-up in `writeup/`.

## Results (real course server)

| | density (bits/base) | accuracy | score |
|---|---|---|---|
| starter baseline (SAMPLE) | 1.44 | 0.92 | 0.87 |
| ours, 10 test passages (250–4000 chars) | 3.49 | 1.000 | **3.49** |
| ours, 104 random excerpts (150–4000 chars), three runs | ~3.4 | 98/104 exact | 3.30, 2.89, 3.26 |

Failures are whole-passage losses when more strands vanish than there are parity strands
(the channel drops ~3% of strands for good). The planner trades this risk against density on
purpose; see "Risk knob" below.

## What each piece is worth (`python3 ablate.py real|sim`)

Real channel, mean score over 150 / 250 / 1000 / 4000-character excerpts (5 each), one piece off at a time:

| variant | mean score | note |
|---|---|---|
| full design | **3.24** | 0/20 failed (3.1% in 160 simulated) |
| no compression | 1.48 | compression is worth 2.2× |
| zlib instead of context mixing | 2.80 | CM ~15% better at every length |
| GG/CC allowed, 2 bits/base | 0.56 | 13/20 failed: strands must be ~40 bases and many never arrive whole |
| no parity | 3.26 | same on short (planner uses none); 20/40 simulated 4000-char passages fail |
| fixed L=200, r=2 | 2.26 | planner matters most for short passages (1.20 vs 2.68) |
| risk cap `MAX_FAIL=0.01` | 2.77 | simulated failures 3.1% → 0.6% |

Full tables: `results/ablation_*_summary.txt`. Explanations for the check-in: `writeup/walkthrough.md`.

## What the channel does (Part 1 of the notebook)

* A strand breaks after **GG or CC** with p ≈ 0.25, elsewhere p ≈ 0.010 (independent).
* Reads shorter than 8 bases are dropped. ~0.6% errors per base (0.46% substitution, 0.11% insertion).
* ~3% of strands are **never read**. Which ones is fixed per sequence; no motif, end, GC or
  length effect was found.

## Design

compress (adaptive context mixing + arithmetic coder, ~3.5–4 bits/char)
→ fixed-length strands `[index | n_data | payload | CRC-8]`
→ each strand unranked into a string with **no GG / CC** (1.83 bits/base)
→ Reed–Solomon (Cauchy, GF(256)) parity strands across strands.
The encoder chooses strand length L and parity count r per passage to maximise
expected score. The decoder finds L from the read-length spike, clusters whole-strand reads,
takes a per-position majority, checks the CRC, rebuilds lost strands, and decompresses.

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
| `results/` | ablation CSVs (one row per run) and summary tables |
| `texts/` | test passages (`python3 texts/make_texts.py` regenerates them) |
| `experiments/` | prototypes (PPM vs context mixing, tuning) and raw channel probes |
| `writeup/` | one-page write-up (HTML source + PDF) and `walkthrough.md` (check-in prep) |

## Workflow

Edit `solution.py`, then `python3 bench.py` and `python3 bench.py real`, then
`python3 make_notebook.py` to refresh the notebook's solution cell (this rebuilds the notebook
without outputs, including Part 4 from `results/ablation_real.csv`). Re-execute with
`jupyter nbconvert --to notebook --execute --inplace ENG_27_Project.ipynb` (needs numpy, matplotlib, requests).

## Risk knob

`MAX_FAIL` in `solution.py` caps the per-passage chance of an unrecoverable loss. At the
default (1.0) the encoder maximises expected score and short passages carry no parity
(~7–10% risk). Capping at 0.01 makes passages near-certain but lowers the expected score,
e.g. for 150-character passages from 2.49 to 1.91.

## Open question for the instructor

We assumed `solution.py` may not ship a pre-trained English model (everything the decoder
needs travels in the strands). If a fixed prior *is* allowed, short passages would compress
much better. Worth asking.
