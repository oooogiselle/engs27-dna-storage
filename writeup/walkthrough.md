# Check-in walkthrough: every decision, why we made it, and where it stops helping

For both of us to read before the check-in. The grade is about explaining engineering decisions,
**when they help and when they don't**, so each section ends with that. Numbers come from
`results/ablation_real_summary.txt` (real server, 5 passages per cell) and
`results/ablation_sim_summary.txt` (simulator, 40 per cell). Rerun with `python3 ablate.py real` / `sim`.

| what we turned off | real mean score | vs full design (3.24) |
|---|---|---|
| compression (raw bytes) | 1.48 | −54% |
| context mixing → zlib | 2.80 | −14% |
| GG/CC-free code | 0.56 | −83% (13 of 20 passages failed) |
| parity strands | 3.26 | same on average, but 20/40 long passages fail in sim |
| per-passage planner → fixed L=200, r=2 | 2.26 | −30% |
| risk capped at 1% | 2.77 | −15%, failures 3.1% → 0.6% (sim) |

---

## The 1-minute version (say this first)

> We measured the channel first. Strands break mostly after GG or CC, about 3% of strands never come back,
> and reads have about 0.6% errors per base. So: we compress the text as hard as we can, because
> with 4 reads per character coverage is plentiful and score is basically density. We write the bits
> with a code that never uses GG or CC, so strands of 200–400 bases survive whole. Each strand carries an
> index and a checksum. We add Reed–Solomon parity strands so lost strands can be rebuilt. And for each
> passage the encoder picks the strand length and number of parity strands that maximise expected score.
> The decoder groups whole-length reads, votes letter by letter, checks the checksum, rebuilds missing
> strands and decompresses. Turning each piece off one at a time shows what it is worth.

---

## 1. What we measured about the channel

**What:** we sent strands we knew, then lined the reads up against them to see exactly where each
molecule broke and what errors it had.

**Findings:**
- **Breaks happen after GG or CC.** The chance of a cut after GG/CC is about 0.25 per site, about 25× the
  0.01 elsewhere. Only the two bases *to the left* matter (we checked the right side and a third base).
- **Elsewhere, breaks are independent.** For strands with no GG/CC, the fraction of reads that are the
  whole strand falls as (1−p)^(L−1) with p ≈ 0.010: 0.75, 0.56, 0.30, 0.13 at L = 30, 60, 120, 200.
  So fragment lengths are geometric once GG/CC is gone.
- **Errors:** ~0.46% substitutions and ~0.11% insertions per base, almost no deletions.
- **~3% of strands are never read at all.** Resent, the same strands vanish again (89 of 90), so it is a
  property of the sequence, but we found no pattern (length, GC content, ends, short motifs), so we treat it
  as random.

**Likely question:** *"How do you know the break depends on two bases and not three?"* We compared the
cut rate grouped by 2 left bases vs 3 left bases. Adding the third base didn't change the rates, and the
right-hand bases didn't matter either (notebook Part 1.2, Fig. 1 in the write-up).

---

## 2. Compression (context mixing + arithmetic coding)

**What it does, plainly:** before writing anything into DNA, we squeeze the text. The compressor predicts
each next bit of the text from what came before (the last 1, 2, 3, 4 or 6 letters, and the current word),
blends those predictions, and spends fewer bits on likely letters. The decoder makes exactly the same
predictions as it goes, so we don't have to send any dictionary or code table. Instead of storing the
length, we code a "more text?" flag before each byte.

**Technical names:** adaptive context-mixing model (orders 0–4, 6 + word hash, logistic mixer, APM) driving a
binary arithmetic coder (Lectures 9–10).

**Why:** with 4 reads per character, there is lots of coverage, so errors are rare and
score ≈ density = 8 × characters / bases. Every bit saved in compression is a direct gain.

**What it's worth:** without compression the mean score falls from 3.24 to 1.48 (2.2× worse). zlib gets
2.80, so context mixing is ~15% better than zlib at every length.

**When it doesn't help / costs:**
- Short passages compress worst, because the model starts knowing nothing about English: density 2.7 at
  150 characters vs 4.1 at 4,000.
- It's slower than zlib (still under a second to decode 4,000 characters).
- We chose not to ship a pre-trained English model (we weren't sure it's allowed). That would help short
  passages most. Worth asking.

**Likely questions:**
- *"Why not Huffman?"* Huffman needs the code table sent, which is expensive for short passages, and it can't
  spend less than 1 bit per letter. Adaptive arithmetic coding sends no table and uses fractional bits.
- *"How does the decoder know where the text ends?"* The "more text?" flag before each byte; it costs a
  fraction of a bit per character and replaces a fixed length header.

---

## 3. The GG/CC-free code (constrained code)

**What it does, plainly:** we never write GG or CC into a strand, because that's where the channel cuts.
We count how many DNA strings of length L have no GG or CC (about 3.56^L of them), number them in order,
and store our bits as "the strand with number X". That gives 1.83 bits per base instead of 2.

**Technical name:** enumerative coding (rank/unrank) of a constrained code (Lecture 18).

**Why:** without GG/CC, breaks drop about 4×, so long strands (200–400 bases) still get many whole reads.

**What it's worth:** this was the most important piece. Without it (plain 2 bits/base), strands have to be
32–52 bases long. The index and checksum then take 23–38% of each strand, and 13 of 20 passages failed
on the real server (mean 0.56).

**Why it failed so badly (good thing to say out loud):** random strands differ. One with several GG/CC
sites almost never arrives whole. In one 1,000-character test, 16 of 78 strands got *zero* whole reads,
more than the 11 parity strands could rebuild. Our simulator used one average break rate for every strand,
so it predicted a score of 2.2; the real score was 0.56. A model is only trustworthy where you've measured it.

**When it doesn't help / costs:** it always costs 8.6% of density per base. It would matter less if we
assembled partial reads (then strands don't need to arrive whole), or for very short strands where few
GG/CC sites occur anyway.

**Likely question:** *"Why not also ban long runs like AAAA?"* We measured: the cut rate depends only on the two
bases to the left, and only GG and CC stand out, so banning anything else would cost density for nothing.

---

## 4. Strand format and decoding

**What:** each strand = `[index (8 bits) | number of data strands (8 bits) | payload | CRC-8]`.
The decoder:
1. finds the strand length L from the spike in the read-length histogram (whole-strand reads all have length L);
2. groups reads of length L that are nearly identical (Hamming distance);
3. takes a letter-by-letter majority vote in each group, which removes the 0.6% errors;
4. converts back to bits and keeps only strands whose CRC checks out;
5. if two candidates claim the same index, keeps the one more reads support.

**Why:** the index lets us put strands in order (reads arrive shuffled). The CRC catches a wrong vote or a
misread strand so it is treated as missing rather than corrupting the text (the parity can then rebuild it).

**When it doesn't help:** most reads are partial and get thrown away. That is what limits how long strands
can be. Assembling overlapping partial reads (de Bruijn graph, Lecture 15) is the next step up, but complex.

**Likely questions:**
- *"Why only 8 bits of CRC?"* It's 8 bits out of ~500 per strand, and a bad strand slipping through by chance
  (1 in 256) also has to beat the correct candidate's vote count. More CRC bits would cost density on every strand.
- *"Why 8-bit index?"* Up to 255 strands covers passages of tens of thousands of characters at our strand
  lengths (we tested 20,000 characters: 59 strands).

---

## 5. Reed–Solomon parity strands

**What it does, plainly:** for n data strands, we add r extra strands computed from all of them, so that
*any* r missing strands can be rebuilt (like RAID for disks).

**Technical name:** systematic Cauchy Reed–Solomon erasure code over GF(256) (Lecture 11).

**Why:** ~3% of strands are lost no matter what, and some get too few reads. A 4,000-character passage
has ~25–30 strands, so losing one is normal.

**What it's worth:** on long passages, without parity, 20 of 40 simulated 4,000-character passages failed
(0 of 40 with parity). Parity costs ~11% of density there, which is a good trade.

**When it doesn't help:**
- **Short passages:** 1–3 strands total. One parity strand adds 33–100% more bases to remove a ~7% risk,
  so the expected score drops. The planner uses no parity for them.
- **~1,000 characters:** break-even. Without parity the score is higher when nothing is lost and zero when
  something is, and on average it comes out the same (simulated 3.17 vs 3.15).

**Likely question:** *"Why not just repeat every strand twice?"* That doubles the bases (halves density) and
still fails if both copies are lost. r parity strands protect against any r losses for r extra strands.

---

## 6. The planner (choosing strand length L and parity count r per passage)

**What it does:** for each passage, the encoder tries every L from 60 to 800 and r from 0 to 12 and picks
the one with the best expected score: density × P(at most r strands fail). A strand "fails" if it's lost
(3.5%) or gets fewer than 8 whole reads (a Poisson count from the coverage, Lectures 6–7).

**Why:** the right L and r depend on passage length. Long passages can afford parity; short ones can't.
Longer strands waste less on index/CRC but get fewer whole reads.

**What it's worth:** a fixed L = 200, r = 2 scores 2.26 instead of 3.24. The gap is biggest on short
passages (1.20 vs 2.68 at 150 characters), because a fixed design adds two parity strands to a passage
that needs one strand.

**When it doesn't help:** at medium lengths a well-chosen fixed design is nearly as good (2.86 vs 3.26 at
1,000 characters). And the planner is only as good as its model: it predicted 5.1% failures where we saw
3.1% in simulation (slightly cautious), but for plain-DNA strands its model was badly wrong (see §3).

---

## 7. The risk knob (`MAX_FAIL`)

**What:** by default we maximise *expected* score, which accepts a few percent chance of losing a passage.
Setting `MAX_FAIL = 0.01` forces the planner to keep failure risk under 1%.

**What it costs:** failures drop from 3.1% to 0.6% (simulated), but the mean score drops from 3.24 to 2.77,
mostly on short passages (2.68 → 1.93 at 150 characters), because they need parity strands they can barely afford.

**When to pick which:** if a single failed passage is a disaster (e.g. only three graded passages), the cap
makes sense; if you care about the average over many passages, expected score wins. Either is defensible;
the point is that it's a choice and we measured the trade.

---

## Things to be honest about

- The simulator (`sim.py`) was fitted to our measurements and matched the full design well, but it was
  far too optimistic for unconstrained strands. We only trust it where we've checked it against the server.
- Real-server runs are 5 per cell, so their failure counts are noisy; failure *rates* come from 40 simulated
  runs per cell.
- We didn't try partial-read assembly or a pre-trained English model; those are the two biggest remaining gains.
