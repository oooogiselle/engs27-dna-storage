# Check-in walkthrough: every decision, why we made it, and where it stops helping

For both of us to read before the check-in. The grade is about explaining engineering decisions,
**when they help and when they don't**, so each section ends with that. Numbers come from
`results/ablation_real_summary.txt` (real server, 5 passages per cell) and
`results/ablation_sim_summary.txt` (simulator, 40 per cell). Rerun with `python3 ablate.py real` / `sim`.

Every technique in `solution.py` is from the course: Huffman coding (3.2.5), parity checks and Hamming
codes (3.3.3–3.3.4), counting writable strings (6.1.4), binomial and Poisson models (2.4, 7.1).

| what we changed | real mean score | vs full design (2.46) |
|---|---|---|
| compression off (raw bytes) | 1.50 | −39% |
| per-passage table → fixed English table | 2.69 | +9% (all of it on short passages) |
| GG/CC-free code off | 0.20 | −92% (17 of 20 passages failed) |
| parity strands off | 2.04 | −17% (4 of 5 4,000-character passages failed) |
| per-passage planner → fixed L=200, r=2 | 1.41 | −43% |
| risk capped at 1% | 1.98 | −20%, simulated failures 5.6% → 0% |
| *reference, not course material:* context mixing instead of Huffman | 3.20 | +30% |

On the ten test passages (`python3 bench.py real`) the full design scored **2.66**, all ten decoded exactly.

---

## The 1-minute version (say this first)

> We measured the channel first. Strands break mostly after GG or CC, about 3% of strands never come back,
> and reads have about 0.6% errors per base. So: we compress the text with a Huffman code built for the
> passage, because with 4 reads per character coverage is plentiful and score is basically density. We write
> the bits with a code that never uses GG or CC, so strands of 250–300 bases survive whole. Each strand carries
> an index and 8 parity bits. We add XOR parity strands, grouped like a Hamming code, so lost strands can be
> rebuilt. And for each passage the encoder picks the strand length and number of parity strands that maximise
> expected score. The decoder groups whole-length reads, votes letter by letter, checks the parity bits,
> rebuilds missing strands and Huffman-decodes. Turning each piece off one at a time shows what it is worth.

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

## 2. Compression (Huffman code built for each passage)

**What it does, plainly:** before writing anything into DNA, we squeeze the text. We count how often each
character appears in *this* passage and build a Huffman code: common characters (space, e, t) get short
codewords, rare ones long codewords. One extra symbol, "end of text", marks where the text stops, so we
don't need a length field.

**How the decoder gets the code:** we send it first, compactly. Only each character's *codeword length* is
sent (4 bits each), plus a 128-bit table saying which characters appear (or a list of characters, if that's
shorter). Both sides then hand out codewords in the same order (shortest first, ties by character), which
rebuilds exactly the same prefix-free code from the lengths alone.

**Technical names:** Huffman coding (Section 3.2.5), sent as a canonical code (lengths only).

**Why:** with 4 reads per character, there is lots of coverage, so errors are rare and
score ≈ density = 8 × characters / bases. Every bit saved in compression is a direct gain.

**What it's worth:** without compression the mean score falls from 2.46 to 1.50 (Huffman is worth 1.6×).
It gets English down to ~4.5 bits/char on long passages, ~6 on 150-character ones, against 8 for raw bytes.

**When it doesn't help / costs:**
- **The table is expensive on short passages.** At 150 characters it costs 1.7 bits per character, 28% of
  everything we store; at 4,000 characters only 0.1. A built-in English table (our fixed-table prototype in
  `partner/solution.py`) sends nothing and scores 2.91 instead of 2.29 at 150 characters, 2.76 vs 2.43 at 250,
  and the same from 1,000 on. The catch: a fixed table is tuned for English prose and would do badly on
  text that isn't (code, other languages, lots of digits); the per-passage table adapts to anything.
- **Huffman codes each character on its own.** It can't use the fact that "q" is followed by "u", or that
  words repeat, and it can't spend less than 1 bit on a character. So it stalls around 4.5 bits/char.
  Methods that use context do better on long passages: zlib reaches 3.7, and the context-mixing coder from
  our earlier version reaches 3.1 bits/char (real score 3.97 vs 2.64 at 4,000 characters). Neither is
  course material, so we kept Huffman.

**Likely questions:**
- *"Why not a fixed English table?"* It's better on short English passages (see above), but the
  per-passage table works for any text and costs almost nothing on long passages. If the graded passages are
  short English, the fixed table would score higher; that's the trade.
- *"How does the decoder know where the text ends?"* The end-of-text symbol: it gets its own codeword,
  which costs a few bits once instead of a fixed length field.
- *"Why 4 bits per codeword length?"* Lengths up to 15 bits. If a code ever came out longer (only for very
  skewed counts), we halve the counts and rebuild until it fits.

---

## 3. The GG/CC-free code (constrained code)

**What it does, plainly:** we never write GG or CC into a strand, because that's where the channel cuts.
We count how many DNA strings of length L have no GG or CC (about 3.56^L of them), number them in order,
and store our bits as "the strand with number X". That gives 1.83 bits per base instead of 2.

**Technical name:** counting writable strings (Section 6.1.4), used as an enumerative code (rank/unrank).

**Why:** without GG/CC, breaks drop about 4×, so long strands (250–300 bases) still get many whole reads.

**What it's worth:** this was the most important piece. Without it (plain 2 bits/base), strands have to be
36–64 bases long. The 28 header and check bits then take 22–39% of each strand, and 17 of 20 passages failed
on the real server (mean 0.20).

**Why it failed so badly (good thing to say out loud):** random strands differ. One with several GG/CC
sites almost never arrives whole, so many strands get *zero* whole reads, more than the parity strands can
rebuild. Our simulator used one average break rate for every strand, so it was far too optimistic here
(32% of simulated passages failed vs 85% real). A model is only trustworthy where you've measured it.

**When it doesn't help / costs:** it always costs 8.6% of density per base. It would matter less if we
assembled partial reads (then strands don't need to arrive whole), or for very short strands where few
GG/CC sites occur anyway.

**Likely question:** *"Why not also ban long runs like AAAA?"* We measured: the cut rate depends only on the two
bases to the left, and only GG and CC stand out, so banning anything else would cost density for nothing.

---

## 4. Strand format and decoding

**What:** each strand = `[index (8 bits) | data strands n (8 bits) | parity strands r (4 bits) | payload | 8 parity bits]`.
Parity bit j is the parity of every 8th bit starting at bit j, so any burst of up to 8 flipped bits changes
at least one of them. The parity bits are also salted with L, so a strand read at the wrong length rarely passes.
The decoder:
1. finds the strand length L from the spike in the read-length histogram (whole-strand reads all have length L);
2. groups reads of length L that are nearly identical (Hamming distance);
3. takes a letter-by-letter majority vote in each group, which removes the 0.6% errors;
4. converts back to bits and keeps only strands whose parity bits check out;
5. if two candidates claim the same index, keeps the one more reads support;
6. reads n and r from the strands themselves (every strand carries them, so it takes a vote).

**Technical name:** parity checks (Section 3.3.3), interleaved 8 ways.

**Why:** the index lets us put strands in order (reads arrive shuffled). The parity bits catch a wrong vote
or a misread strand so it is treated as missing rather than corrupting the text (the XOR strands can then
rebuild it). Carrying n and r on every strand means nothing is shared outside the DNA.

**When it doesn't help:** most reads are partial and get thrown away. That is what limits how long strands
can be. Assembling overlapping partial reads is the next step up, but complex.

**Likely questions:**
- *"Is 8 parity bits enough?"* It's 8 bits out of ~500 per strand. A wrong strand slips through by chance
  1 time in 256, and it also has to out-vote the correct candidate. More check bits would cost density on
  every strand.
- *"Why not a CRC?"* A CRC catches more error patterns, but it isn't course material; the interleaved
  parity catches every burst of up to 8 flips, which is what a bad vote or a shifted read looks like.

---

## 5. XOR parity strands (grouped like a Hamming code)

**What it does, plainly:** a parity strand is the XOR of a group of data strands. If exactly one strand of
the group is lost, XOR-ing the parity strand with the ones that did arrive gives it back. With r parity
strands, every data strand gets an r-bit label with at least two 1s (Hamming skips the powers of two for the
same reason), and it joins group j when bit j of its label is 1. So every data strand is covered by at least
two parity strands. The decoder repeats one step: find a group missing exactly one strand, rebuild it.
Rebuilding one can leave another group missing only one, and so on.

**Technical name:** parity checks across strands with Hamming-code groups (Section 3.3.4), decoded by
repeated single-erasure repair ("peeling").

**Why:** ~3% of strands are lost no matter what, and some get too few reads. A 4,000-character passage
has ~40 data strands, so losing one or two is normal.

**What it's worth:** on long passages, without parity, 4 of 5 real and 28 of 40 simulated 4,000-character
passages failed; with it, 0 of 5 real and 2 of 40 simulated.

**When it doesn't help / costs:**
- **It needs more parity strands than the best possible code.** Two lost strands that sit in the same
  groups can't be untangled. So at 4,000 characters the planner uses 6–7 parity strands on 36–58 data strands
  (density 3.07 → 2.64, −14%). A Reed–Solomon code, which rebuilds *any* r lost strands, needed only 2–3 in
  our earlier version. That's a real cost of keeping to course methods.
- **Short passages:** 1–3 strands total. One parity strand adds 33–100% more bases to remove a ~10% risk,
  so the expected score drops. The planner uses no parity for them.
- **~1,000 characters:** close to break-even with one parity strand (simulated 2.62 with vs 2.45 without,
  where 7 of 40 failed; on the real server none of 5 lost a strand, so no-parity happened to score higher).
- **r = 2 is useless**: the only 2-bit label with two 1s is 11, so both parity strands would cover the same
  strands. The planner skips it (that's also why the fixed L = 200, r = 2 design failed 3 of 5 long passages).

**Likely question:** *"Why not just repeat every strand twice?"* That doubles the bases (halves density) and
still fails if both copies are lost. 6–7 parity strands on 40 data strands is ~16% extra.

---

## 6. The planner (choosing strand length L and parity count r per passage)

**What it does:** for each passage, the encoder tries every L from 60 to 800 and r from 0 to 12 and picks
the one with the best expected score: density × P(every data strand is recovered).

**How P is worked out:** a strand "fails" if it's lost (3.5%) or gets fewer than 8 whole reads (a Poisson
count from the coverage, Section 7.1). The number of failed strands out of n + r is binomial (Section 2.4).
For each number k of failed strands, we work out the fraction of ways to lose k strands that the XOR groups
can still repair (every pattern for k ≤ 2 when there are at most 800 of them, otherwise 120 random
patterns with a fixed seed, so the encoder's choice is repeatable), and add up
P(k failed) × P(repairable | k failed) (law of total probability).

**Why:** the right L and r depend on passage length. Long passages can afford parity; short ones can't.
Longer strands waste less on the header but get fewer whole reads.

**What it's worth:** a fixed L = 200, r = 2 scores 1.41 instead of 2.46. The gap is biggest on short
passages (1.20 vs 2.29 at 150 characters), because a fixed design adds two parity strands to a passage
that needs one strand, and on long ones, where r = 2 is too little (3 of 5 failed).

**When it doesn't help:** the planner is only as good as its model. For our design it is close (predicted
6.5% failures, 5.6% seen in 160 simulated passages), but for plain-DNA strands its model was badly wrong (see §3).

---

## 7. The risk knob (`MAX_FAIL`)

**What:** by default we maximise *expected* score, which accepts a few percent chance of losing a passage.
Setting `MAX_FAIL = 0.01` forces the planner to keep failure risk under 1%.

**What it costs:** simulated failures drop from 5.6% to 0%, but the mean score drops from 2.46 to 1.98,
mostly on short passages (2.29 → 1.65 at 150 characters), because they need parity strands they can barely afford.

**When to pick which:** if a single failed passage is a disaster (e.g. only three graded passages), the cap
makes sense; if you care about the average over many passages, expected score wins. Either is defensible;
the point is that it's a choice and we measured the trade.

---

## Things to be honest about

- **We gave up score to stay within course methods.** Our earlier version (context mixing + Reed–Solomon +
  CRC-8) scored 3.24 on the same real-channel ablation; this one scores 2.46. The ablation's context-mixing
  row (3.20) shows almost all of that gap is the compressor; the rest is XOR needing more parity strands.
  The earlier version is kept in `experiments/solution_cm.py`.
- The simulator (`sim.py`) was fitted to our measurements and matched the full design well, but it was
  far too optimistic for unconstrained strands. We only trust it where we've checked it against the server.
- Real-server runs are 5 per cell, so their failure counts are noisy; failure *rates* come from 40 simulated
  runs per cell.
- Short passages carry a ~10% risk of total loss by design (no parity): 5 of 40 simulated 150-character
  passages failed, and the notebook's simulator run on SAMPLE (4 strands) lost one and scored 0.
  On the real server 53 of 54 runs were exact; the one failure was a 1,000-character excerpt in the
  notebook's stress test, where more strands were lost than its parity groups could repair.
- We didn't try partial-read assembly; it's the biggest remaining gain within the course material, along
  with a fixed English Huffman table for short passages.
