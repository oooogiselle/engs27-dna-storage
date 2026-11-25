"""ENGS 27 final project: storing English text in DNA.

Pipeline
  encode: text -> Huffman code built for this passage (its code lengths are sent first)
          -> split into strands of equal length L, each strand = [index | n | r | payload | parity]
          -> each strand's bits are mapped to a DNA string with NO "GG" and NO "CC"
             (enumerative code; the real channel breaks ~25x more often after GG/CC)
          -> r extra XOR parity strands, grouped like a Hamming code, rebuild lost strands
  decode: estimate L from the spike in the read-length histogram (full-length reads)
          -> cluster full-length reads by Hamming distance, per-position majority vote
          -> unrank each consensus, check its parity bits, keep the best-supported copy
          -> rebuild lost strands from the XOR parity strands
          -> concatenate payloads and Huffman-decode; stop at the first strand still missing.

Everything the decoder needs travels in the strands: the Huffman code lengths are the
first bits of the message, and every strand carries n (data strands) and r (parity strands).
Every technique here is from the course: Huffman coding (3.2.5), parity checks and Hamming
codes (3.3.3-3.3.4), counting writable strings (6.1.4), binomial and Poisson models (2.4, 7.1).
"""
import collections
import heapq
import math
import random

# ----------------------------------------------------------------------------------
# Tunable parameters
# ----------------------------------------------------------------------------------
CHECK_BITS = 8     # parity bits at the end of every strand
HDR_BITS = 20      # 8-bit strand index + 8-bit number of data strands + 4-bit parity count
P_LOST = 0.035     # measured ~3% of strands never come back; plan for a bit more
P_BREAK = 0.0101   # measured break probability per site once GG/CC are banned
L_RANGE = range(60, 801, 4)
MIN_FULL = 8       # want at least this many full-length reads per strand
MAX_FAIL = 1.0     # cap on P(passage unrecoverable); 1.0 = pure expected-score maximising
MAX_R = 12         # most parity strands (fits in the 4-bit field)

# ----------------------------------------------------------------------------------
# 1. Compressor: a Huffman code built for this passage (Section 3.2.5)
# ----------------------------------------------------------------------------------
# Symbols are the passage's bytes plus an end-of-text symbol (EOT), so no length field is
# needed. The decoder needs the code, so we send it first, compactly: each symbol's codeword
# length. Both sides then hand out codewords in the same order (shortest first, then by
# symbol), which gives the same prefix-free code from the lengths alone.
EOT = 256
MAX_LEN = 15       # codeword lengths are sent in 4 bits


def _huffman_lengths(freq):
    """Codeword length of every symbol: repeatedly merge the two rarest groups; every
    merge makes each symbol in them one bit longer."""
    heap = [(f, s, (s,)) for s, f in freq.items()]
    heapq.heapify(heap)
    if len(heap) == 1:
        return {heap[0][1]: 1}
    depth = collections.Counter()
    tie = 1000
    while len(heap) > 1:
        f1, _, a = heapq.heappop(heap)
        f2, _, b = heapq.heappop(heap)
        for s in a + b:
            depth[s] += 1
        heapq.heappush(heap, (f1 + f2, tie, a + b))
        tie += 1
    return dict(depth)


def _code_lengths(data):
    freq = collections.Counter(data)
    freq[EOT] = 1
    lens = _huffman_lengths(freq)
    while max(lens.values()) > MAX_LEN:           # very rare: flatten the counts and retry
        freq = {s: (f + 1) // 2 for s, f in freq.items()}
        lens = _huffman_lengths(freq)
    return lens


def _codebook(lens):
    """Codewords from lengths: shortest first, ties by symbol (the same on both sides)."""
    code, prev_len, book = 0, 0, {}
    for s in sorted(lens, key=lambda s: (lens[s], s)):
        code <<= lens[s] - prev_len
        book[s] = (code, lens[s])
        code += 1
        prev_len = lens[s]
    return book


def _write_header(lens):
    """Send the code lengths. Two layouts; a flag bit says which one is shorter:
    0: 9-bit symbol count, then (8-bit symbol, 4-bit length) per byte symbol;
    1: 128-bit table of which ASCII characters appear, then a 4-bit length for each.
    Both end with the 4-bit length of EOT."""
    syms = sorted(s for s in lens if s != EOT)
    a = [0] + _to_bits(len(syms), 9)
    for s in syms:
        a += _to_bits(s, 8) + _to_bits(lens[s], 4)
    best = a
    if all(s < 128 for s in syms):
        b = [1] + [int(c in lens) for c in range(128)]
        for s in syms:
            b += _to_bits(lens[s], 4)
        if len(b) < len(a):
            best = b
    return best + _to_bits(lens[EOT], 4)


def _read_header(bits):
    """Inverse of _write_header. Returns (lens, position after the header)."""
    pos, lens = 1, {}
    if bits[0] == 0:
        count = _to_int(bits[pos:pos + 9])
        pos += 9
        for _ in range(count):
            lens[_to_int(bits[pos:pos + 8])] = _to_int(bits[pos + 8:pos + 12])
            pos += 12
    else:
        present = [c for c in range(128) if bits[pos + c]]
        pos += 128
        for c in present:
            lens[c] = _to_int(bits[pos:pos + 4])
            pos += 4
    lens[EOT] = _to_int(bits[pos:pos + 4])
    pos += 4
    if any(not 1 <= v <= MAX_LEN for v in lens.values()):
        raise ValueError("damaged header")
    return lens, pos


def compress(data):
    lens = _code_lengths(data)
    book = _codebook(lens)
    bits = _write_header(lens)
    for s in list(data) + [EOT]:
        code, n = book[s]
        bits += _to_bits(code, n)
    return bits


def decompress(bits, stop_bit=None, limit=1 << 20):
    """Decode bytes. If stop_bit is given, stop before reading past that bit position
    (the data there is missing, so later output would be garbage)."""
    lens, pos = _read_header(bits)
    lookup = {(n, code): s for s, (code, n) in _codebook(lens).items()}
    end = len(bits) if stop_bit is None else min(len(bits), stop_bit)
    out = []
    while len(out) < limit:
        code, n = 0, 0
        while True:
            if pos >= end or n >= MAX_LEN:
                return bytes(out)
            code = code * 2 + bits[pos]
            pos += 1
            n += 1
            s = lookup.get((n, code))
            if s is not None:
                break
        if s == EOT:
            break
        out.append(s)
    return bytes(out)


# ----------------------------------------------------------------------------------
# 2. Constrained code: DNA strings with no GG and no CC, ranked as integers
# ----------------------------------------------------------------------------------
BASES = "ACGT"
_COUNT_CACHE = {}


def _ok(prev, b):
    return not (prev == b and b in "GC")


def _counts(n):
    """N[k][prev] = number of valid strings of length k that may follow base prev
    (prev = None at the start)."""
    if n not in _COUNT_CACHE:
        N = [{p: 1 for p in BASES + "S"}]
        for k in range(1, n + 1):
            N.append({p: sum(N[k - 1][b] for b in BASES if p == "S" or _ok(p, b))
                      for p in BASES + "S"})
        _COUNT_CACHE[n] = N
    return _COUNT_CACHE[n]


def capacity(L):
    """Bits one constrained strand of length L can carry."""
    return _counts(L)[L]["S"].bit_length() - 1


def unrank(x, L):
    N, prev, out = _counts(L), "S", []
    for k in range(L, 0, -1):
        for b in BASES:
            if prev != "S" and not _ok(prev, b):
                continue
            c = N[k - 1][b]
            if x < c:
                out.append(b)
                prev = b
                break
            x -= c
    return "".join(out)


def rank(s):
    L = len(s)
    N, prev, x = _counts(L), "S", 0
    for i, ch in enumerate(s):
        if ch not in BASES or (prev != "S" and not _ok(prev, ch)):
            return None
        for b in BASES:
            if b == ch:
                break
            if prev == "S" or _ok(prev, b):
                x += N[L - i - 1][b]
        prev = ch
    return x


# ----------------------------------------------------------------------------------
# 3. XOR parity strands across strands, grouped like a Hamming code (Sections 3.3.3-3.3.4)
# ----------------------------------------------------------------------------------
# A parity strand is the XOR of a group of data strands: if exactly one strand of the group
# is lost, XOR-ing the parity strand with the strands that did arrive gives it back.
# With r parity strands, every data strand gets a label: an r-bit number with at least two
# 1s (Hamming skips the powers of two for the same reason). Data strand i is in group j when
# bit j of its label is 1, so each data strand is covered by at least two parity strands.
# Labels with two 1s are used first, then three, ...
# Decoding repeats one step: find a group missing exactly one strand and rebuild it.
def _labels(r):
    if r == 1:
        return [1]                       # one parity strand: the XOR of all data strands
    return sorted((x for x in range(1, 2 ** r) if bin(x).count("1") >= 2),
                  key=lambda x: (bin(x).count("1"), x))


def _groups(n, r):
    """groups[j] = the data strands that parity strand j covers."""
    if r == 0:
        return []
    labels = _labels(r)
    groups = [[] for _ in range(r)]
    for i in range(n):
        label = labels[i % len(labels)]
        for j in range(r):
            if label >> j & 1:
                groups[j].append(i)
    return groups


def _xor(rows, K):
    out = [0] * K
    for row in rows:
        for t in range(K):
            out[t] ^= row[t]
    return out


def _parity(data, r):
    K = len(data[0])
    return [_xor([data[i] for i in g], K) for g in _groups(len(data), r)]


def _recover(data, parity, n, r):
    """data: {i: payload} for received data strands; parity: {j: payload}.
    Rebuilds what it can and returns the updated {i: payload} (may still have gaps)."""
    data = dict(data)
    groups = _groups(n, r)
    progress = True
    while progress and len(data) < n:
        progress = False
        for j, g in enumerate(groups):
            if j not in parity:
                continue
            missing = [i for i in g if i not in data]
            if len(missing) == 1:
                K = len(parity[j])
                data[missing[0]] = _xor([parity[j]] + [data[i] for i in g if i != missing[0]], K)
                progress = True
    return data


# How likely is the passage to be rebuilt? For n data and r parity strands, _fixable(n, r, k)
# is the fraction of ways to lose k of the N = n + r strands that the decoder survives
# (exact for k <= 2, Monte Carlo with a fixed seed for k >= 3). Combined with the binomial
# chance of losing k strands (law of total probability), it gives P(success).
_FIX_CACHE = {}


def _survives(n, masks, lost):
    """Peeling on bitmasks: data strands are bits 0..n-1, parity strand j is bit n+j."""
    missing = lost & ((1 << n) - 1)
    live = [g for j, g in enumerate(masks) if not lost >> (n + j) & 1]
    progress = True
    while missing and progress:
        progress = False
        for g in live:
            x = g & missing
            if x and x & (x - 1) == 0:   # exactly one strand of this group is missing
                missing ^= x
                progress = True
    return missing == 0


def _fixable(n, r, k):
    key = (n, r, k)
    if key not in _FIX_CACHE:
        N = n + r
        masks = [sum(1 << i for i in g) for g in _groups(n, r)]
        if k == 0:
            f = 1.0
        elif r == 0:
            f = 0.0
        elif k == 1:
            f = sum(_survives(n, masks, 1 << a) for a in range(N)) / N
        elif k == 2 and N * (N - 1) // 2 <= 800:
            f = sum(_survives(n, masks, (1 << a) | (1 << b))
                    for a in range(N) for b in range(a + 1, N)) / (N * (N - 1) // 2)
        else:
            rng = random.Random(7919 * N + 31 * r + k)
            trials, ok = 120, 0
            for _ in range(trials):
                lost = 0
                while bin(lost).count("1") < k:
                    lost |= 1 << int(rng.random() * N)
                ok += _survives(n, masks, lost)
            f = ok / trials
        _FIX_CACHE[key] = f
    return _FIX_CACHE[key]


def _p_success(n, r, q):
    """P(every data strand is recovered) when each of the n + r strands fails with prob q."""
    N = n + r
    total, mass, term = 0.0, 0.0, (1 - q) ** N
    for k in range(N + 1):
        total += term * _fixable(n, r, k)
        mass += term
        if 1 - mass < 1e-5 or k >= 10:
            break
        term = term * (N - k) / (k + 1) * q / (1 - q) if q < 1 else 0.0
    return min(total, 1.0)


# ----------------------------------------------------------------------------------
# 4. Strand framing: [index 8 | n_data 8 | r 4 | K payload bytes | pad | 8 parity bits]
# ----------------------------------------------------------------------------------
def _check(bits, salt=0):
    """8 parity bits (Section 3.3.3): bit j is the parity of every 8th bit starting at j.
    Any burst of up to 8 flipped bits changes at least one of them. The salt (from L)
    makes a strand read at the wrong length unlikely to pass."""
    p = [(salt >> (CHECK_BITS - 1 - j)) & 1 for j in range(CHECK_BITS)]
    for i, b in enumerate(bits):
        p[i % CHECK_BITS] ^= b
    return p


def _to_int(bits):
    x = 0
    for b in bits:
        x = x * 2 + b
    return x


def _to_bits(x, n):
    return [(x >> (n - 1 - i)) & 1 for i in range(n)]


def _K(L):
    """Payload bytes per strand of length L."""
    return (capacity(L) - HDR_BITS - CHECK_BITS) // 8


def _frame(index, n, r, payload, L):
    cap = capacity(L)
    bits = _to_bits(index, 8) + _to_bits(n, 8) + _to_bits(r, 4)
    for byte in payload:
        bits += _to_bits(byte, 8)
    bits += [0] * (cap - CHECK_BITS - len(bits))
    bits += _check(bits, L & 0xFF)
    return unrank(_to_int(bits), L)


def _unframe(strand):
    L = len(strand)
    x = rank(strand)
    if x is None:
        return None
    cap = capacity(L)
    if x >> cap:
        return None
    bits = _to_bits(x, cap)
    if _check(bits[:-CHECK_BITS], L & 0xFF) != bits[-CHECK_BITS:]:
        return None
    index, n, r = _to_int(bits[:8]), _to_int(bits[8:16]), _to_int(bits[16:20])
    K = _K(L)
    payload = [_to_int(bits[HDR_BITS + 8 * t:HDR_BITS + 8 * t + 8]) for t in range(K)]
    return index, n, r, payload


# ----------------------------------------------------------------------------------
# 5. Choosing strand length L and number of parity strands r for each passage
# ----------------------------------------------------------------------------------
def _plan(nbits, nchars):
    best = None
    for L in L_RANGE:
        K = _K(L)
        if K < 4:
            continue
        n = max(1, -(-nbits // (8 * K)))
        if n > 200:
            continue
        for r in [0, 1] + list(range(3, MAX_R + 1)):
            N = n + r
            # reads per surviving strand, and expected full-length reads among them
            lam = 4 * nchars / (N * (1 - P_LOST)) * (1 - P_BREAK) ** (L - 1)
            # P(fewer than MIN_FULL full reads): too few to find L and vote reliably
            term, p_few = math.exp(-lam), 0.0
            for k in range(MIN_FULL):
                p_few += term
                term *= lam / (k + 1)
            q = P_LOST + (1 - P_LOST) * p_few
            ok = _p_success(n, r, q)
            if 1 - ok > MAX_FAIL:
                continue
            ev = ok / (N * L)
            if best is None or ev > best[0]:
                best = (ev, L, n, r)
    if best is None:
        return _plan_fallback(nbits)
    return best[1:]


def _plan_fallback(nbits):
    L = L_RANGE[-1]
    return L, max(1, -(-nbits // (8 * _K(L)))), 0


# ----------------------------------------------------------------------------------
# 6. encode / decode
# ----------------------------------------------------------------------------------
def encode(text):
    bits = compress(text.encode("utf-8"))
    L, n, r = _plan(len(bits), max(1, len(text)))
    K = _K(L)
    bits = bits + [0] * (n * K * 8 - len(bits))
    data = [[_to_int(bits[8 * (i * K + t):8 * (i * K + t) + 8]) for t in range(K)] for i in range(n)]
    strands = [_frame(i, n, r, d, L) for i, d in enumerate(data)]
    strands += [_frame(n + j, n, r, p, L) for j, p in enumerate(_parity(data, r))]
    return strands


def _estimate_L(reads):
    cnt = collections.Counter(len(r) for r in reads)
    best, score = None, -1
    for l, c in cnt.items():
        if l < 12:
            continue
        s = c - cnt.get(l + 2, 0) - cnt.get(l + 3, 0)    # a spike above its neighbours
        if s > score:
            best, score = l, s
    return best


def _consensus_candidates(reads, L):
    """Cluster reads of length exactly L by Hamming distance; return (consensus, support)."""
    full = collections.Counter(r for r in reads if len(r) == L)
    k = 4                                   # pigeonhole: <=3 errors => one chunk exact
    bounds = [(L * j // k, L * (j + 1) // k) for j in range(k)]
    tol = max(3, L // 10)
    clusters = []                           # [rep, Counter-of-reads, total]
    chunk_idx = collections.defaultdict(list)
    for r, c in full.most_common():
        hit, seen = None, set()
        for j, (a, b) in enumerate(bounds):
            for ci in chunk_idx.get((j, r[a:b]), ()):
                if ci in seen:
                    continue
                seen.add(ci)
                if sum(x != y for x, y in zip(clusters[ci][0], r)) <= tol:
                    hit = ci
                    break
            if hit is not None:
                break
        if hit is None:
            hit = len(clusters)
            clusters.append([r, collections.Counter(), 0])
            for j, (a, b) in enumerate(bounds):
                chunk_idx[(j, r[a:b])].append(hit)
        clusters[hit][1][r] += c
        clusters[hit][2] += c
    out = []
    for rep, members, total in clusters:
        if total > 1:
            cols = [collections.Counter() for _ in range(L)]
            for r, c in members.items():
                for i, ch in enumerate(r):
                    cols[i][ch] += c
            cons = "".join(col.most_common(1)[0][0] for col in cols)
            out.append((cons, total))
            if cons != rep:
                out.append((rep, members[rep]))
        else:
            out.append((rep, 1))
    return out


def decode(reads):
    try:
        return _decode(reads)
    except Exception:
        return ""


def _decode(reads):
    reads = [r for r in reads if r]
    if not reads:
        return ""
    # Try the likely strand lengths; keep the one whose reads decode to the most
    # parity-valid strands (weighted by how many reads support them).
    choice = None
    for L in _candidate_lengths(reads):
        best, nvotes = {}, collections.Counter()
        for cons, support in _consensus_candidates(reads, L):
            got = _unframe(cons)
            if got is None:
                continue
            i, n, r, payload = got
            nvotes[(n, r)] += support
            if i not in best or support > best[i][0]:
                best[i] = (support, payload)
        total = sum(v[0] for v in best.values())
        if best and (choice is None or total > choice[0]):
            choice = (total, L, best, nvotes)
    if choice is None:
        return ""
    _, L, best, nvotes = choice
    n, r = nvotes.most_common(1)[0][0]      # every strand carries n and r, so take a vote
    data = {i: v[1] for i, v in best.items() if i < n}
    parity = {i - n: v[1] for i, v in best.items() if n <= i < n + r}
    data = _recover(data, parity, n, r)
    K = _K(L)
    full, stop = [], None
    for i in range(n):                      # decode up to the first strand still missing
        if i not in data:
            stop = 8 * K * i
            break
        full.append(data[i])
    bits = []
    for d in full:
        for byte in d:
            bits += _to_bits(byte, 8)
    out = decompress(bits, stop_bit=stop, limit=max(20000, 100 * len(bits)))
    return out.decode("utf-8", errors="replace")


def _candidate_lengths(reads):
    first = _estimate_L(reads)
    cnt = collections.Counter(len(r) for r in reads if len(r) >= 20)
    out = [first, first - 1, first + 1] + [l for l, _ in cnt.most_common(8)]
    # long lengths that occur several times are strong candidates too (full reads)
    out += [l for l, c in cnt.items() if c >= 3 and l >= 60]
    seen, res = set(), []
    for l in out:
        if l and l >= 20 and l not in seen:
            seen.add(l)
            res.append(l)
    return res
