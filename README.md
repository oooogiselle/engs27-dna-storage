# ENGS 27 final project: storing text in DNA

Hand-in: **`solution.py`** (`encode(text) -> list[str]`, `decode(reads) -> str`) and the
one-page write-up in `writeup/`.

## Results (real course server)

| | density (bits/base) | accuracy | score |
|---|---|---|---|
| starter baseline (SAMPLE) | 1.44 | 0.92 | 0.87 |
| ours, 10 test passages (150–4000 chars) | 3.49 | 1.000 | **3.49** |
| ours, 80 random excerpts over two runs | | 75/80 perfect | ~3.1 mean |

Failures are whole-passage losses when more strands vanish than there are parity strands
(the channel drops ~3% of strands for good). The planner trades this risk against density on
purpose; see "Risk knob" below.

## What the channel does (Part 1 of the notebook)

* A strand breaks after **GG or CC** with p ≈ 0.25, elsewhere p ≈ 0.010 (independent).
* Reads shorter than 8 bases are dropped. ~0.5% errors per base (85% substitution, 15% insertion).
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
| `texts/` | test passages (`python3 texts/make_texts.py` regenerates them) |
| `experiments/` | prototypes (PPM vs context mixing, tuning) and raw channel probes |
| `writeup/` | one-page write-up (HTML source + PDF) and figures |

## Workflow

Edit `solution.py`, then `python3 bench.py` and `python3 bench.py real`, then
`python3 make_notebook.py` to refresh the notebook's solution cell. Re-execute with
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
