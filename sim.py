"""Channel simulator fitted to measurements of the real course server (2026-11-16).

Measured (see exploration/ and the notebook, Part 1):
  * a molecule breaks after a GG or CC dinucleotide with p ~ 0.25, elsewhere p ~ 0.010,
    independently; the bases to the right of the cut do not matter
  * pieces shorter than 8 bases are never returned; a read is a random kept piece
  * ~0.5% errors per base: ~85% substitutions, ~15% insertions, deletions rare
  * ~3% of strands are never read at all. Which ones is fixed per sequence (resending
    gives the same losses) but shows no motif, end, length or GC dependence, so we
    model it as a pseudo-random hash of the sequence (salted by seed).
"""
import random
import zlib

P_HOT, P_COLD, MIN_LEN = 0.25, 0.0101, 8
P_SUB, P_INS, P_DEL = 0.0043, 0.0008, 0.0001
P_LOST = 0.03


class SimChannel:
    def __init__(self, seed=None, p_hot=P_HOT, p_cold=P_COLD, p_sub=P_SUB, p_ins=P_INS, p_del=P_DEL,
                 p_lost=P_LOST):
        self.rng = random.Random(seed)
        self.salt = (seed or 0) * 7919
        self.p_lost = p_lost
        self.p_hot, self.p_cold = p_hot, p_cold
        self.p_sub, self.p_ins, self.p_del = p_sub, p_ins, p_del

    def _one(self, strands):
        rng = self.rng
        while True:
            s = rng.choice(strands)
            cuts = [i for i in range(1, len(s))
                    if rng.random() < (self.p_hot if i >= 2 and s[i - 2] == s[i - 1] and s[i - 1] in "GC" else self.p_cold)]
            pieces = [s[a:b] for a, b in zip([0] + cuts, cuts + [len(s)])]
            pieces = [x for x in pieces if len(x) >= MIN_LEN]
            if pieces:
                return self._noise(rng.choice(pieces))

    def _noise(self, r):
        rng, out = self.rng, []
        for c in r:
            u = rng.random()
            if u < self.p_del:
                continue
            if u < self.p_del + self.p_sub:
                out.append(rng.choice([b for b in "ACGT" if b != c]))
            else:
                out.append(c)
            if rng.random() < self.p_ins:
                out.append(rng.choice("ACGT"))
        return "".join(out)

    def _lost(self, s):
        return (zlib.crc32(s.encode()) ^ self.salt) % 100000 < self.p_lost * 100000

    def reads(self, strands, n):
        strands = [s for s in strands if not self._lost(s)]
        if not strands:
            return []
        return [self._one(strands) for _ in range(n)]
