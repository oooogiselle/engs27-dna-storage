"""ENGS 27 final project: storing English text in DNA.

Pipeline
  encode: text -> context-mixing compressor -> bit stream
          -> split into strands of equal length L, each strand = [index | payload | CRC-8]
          -> each strand's bits are mapped to a DNA string with NO "GG" and NO "CC"
             (enumerative code; the real channel breaks ~25x more often after GG/CC)
  decode: estimate L from the spike in the read-length histogram (full-length reads)
          -> cluster full-length reads by Hamming distance, per-position majority vote
          -> unrank each consensus, check CRC, keep the best-supported copy of each index
          -> concatenate payloads and decompress; stop at the first missing strand.

Nothing about the text is shared outside the strands: the compressor is adaptive and
the decoder rebuilds the same model as it goes.
"""
import collections
import math

# ----------------------------------------------------------------------------------
# Tunable parameters
# ----------------------------------------------------------------------------------
CRC_BITS = 8
HDR_BITS = 16      # 8-bit strand index + 8-bit number of data strands
P_LOST = 0.035     # measured ~3% of strands never come back; plan for a bit more
P_BREAK = 0.0101   # measured break probability per site once GG/CC are banned
L_RANGE = range(60, 801, 4)
MIN_FULL = 8       # want at least this many full-length reads per strand
MAX_FAIL = 1.0     # cap on P(passage unrecoverable); 1.0 = pure expected-score maximising

# ----------------------------------------------------------------------------------
# 1. Compressor: tiny bitwise context mixing (orders 0-4, 6 and a word model)
# ----------------------------------------------------------------------------------
_PREC = 1 << 16


def _stretch(p):
    return math.log(p / (1 - p))


def _squash(x):
    if x > 30:
        x = 30
    elif x < -30:
        x = -30
    return 1 / (1 + math.exp(-x))


class _BitCoder:
    """Binary arithmetic coder, 32-bit range. Encoder when bits is None."""

    def __init__(self, bits=None):
        self.lo, self.hi = 0, 0xFFFFFFFF
        self.bits = bits
        if bits is None:
            self.out = []
        else:
            self.pos, self.x = 0, 0
            for _ in range(32):
                self.x = (self.x << 1) | self._nb()

    def _nb(self):
        b = self.bits[self.pos] if self.pos < len(self.bits) else 0
        self.pos += 1
        return b

    def code(self, p1, bit=None):
        q = min(max(int(p1 * _PREC), 1), _PREC - 1)
        r = self.hi - self.lo
        mid = self.lo + (r >> 16) * q + (((r & 0xFFFF) * q) >> 16)
        if bit is None:
            bit = 1 if self.x <= mid else 0
        if bit:
            self.hi = mid
        else:
            self.lo = mid + 1
        while (self.lo ^ self.hi) & 0x80000000 == 0:
            if self.bits is None:
                self.out.append(self.hi >> 31)
            else:
                self.x = ((self.x << 1) & 0xFFFFFFFF) | self._nb()
            self.lo = (self.lo << 1) & 0xFFFFFFFF
            self.hi = ((self.hi << 1) & 0xFFFFFFFF) | 1
        return bit

    def finish(self):
        for i in range(32):
            self.out.append((self.lo >> (31 - i)) & 1)
            pref = self.lo >> (31 - i) << (31 - i)
            if pref >= self.lo and (pref | ((1 << (31 - i)) - 1)) <= self.hi:
                break
        return self.out


class _Model:
    ORDERS = (1, 2, 3, 4, 6)
    LR0, LR1, LIMIT = 0.01, 0.3, 255     # mixer learning rate decays as the text goes on

    def __init__(self):
        self.t = {}
        self.w = {}
        self.apm = {}
        self.hist = []
        self.wordh = 0
        self.more = 0      # bytes seen, drives the end-of-text flag

    def p_more(self):
        return (self.more + 1) / (self.more + 2)

    def code_byte(self, coder, byte=None):
        h = self.hist
        cs = [("o0", 0)]
        for k in self.ORDERS:
            cs.append((k, bytes(h[-k:])) if len(h) >= k else (k, None))
        cs.append(("w", self.wordh))
        node = 1
        last = h[-1] if h else 0
        lr = self.LR0 + self.LR1 / math.sqrt(1 + len(h))
        for i in range(8):
            bit_in = None if byte is None else (byte >> (7 - i)) & 1
            slots, st = [], []
            for m, c in cs:
                if c is None:
                    slots.append(None)
                    st.append(0.0)
                    continue
                key = (m, c, node)
                e = self.t.get(key)
                if e is None:
                    e = self.t[key] = [0.5, 0]
                slots.append(e)
                st.append(_stretch(min(max(e[0], 1e-4), 1 - 1e-4)) if e[1] else 0.0)
            st.append(0.3)
            sel = (node, min(len(h), 2))
            w = self.w.get(sel)
            if w is None:
                w = self.w[sel] = [0.3] * len(st)
            dot = sum(a * b for a, b in zip(w, st))
            p = _squash(dot)
            bucket = 0 if dot <= -8 else (32 if dot >= 8 else int((dot + 8) * 2))
            ak = (last, node, bucket)
            a = self.apm.get(ak)
            pf = p if a is None else 0.5 * (p + a)
            bit = coder.code(pf, bit_in)
            a0 = p if a is None else a
            self.apm[ak] = a0 + 0.05 * (bit - a0)
            err = bit - p
            for j in range(len(w)):
                w[j] += lr * err * st[j]
            for e in slots:
                if e is not None:
                    e[1] = min(e[1] + 1, self.LIMIT)
                    e[0] += (bit - e[0]) / (e[1] + 0.5)
            node = node * 2 + bit
        b = node & 0xFF
        h.append(b)
        if 65 <= b <= 90 or 97 <= b <= 122 or b >= 128:
            self.wordh = (self.wordh * 0x2F0B3 + (b | 32) + 1) & 0xFFFFFFFFFF   # explicit, version-proof
        else:
            self.wordh = 0
        return b


def compress(data):
    m, c = _Model(), _BitCoder()
    for b in data:
        c.code(m.p_more(), 1)
        m.more += 1
        m.code_byte(c, b)
    c.code(m.p_more(), 0)
    return c.finish()


def decompress(bits, stop_bit=None, limit=1 << 20):
    """Decode bytes. If stop_bit is given, stop once the decoder has had to read past
    that bit position (the data there is missing, so later output would be garbage)."""
    m, c = _Model(), _BitCoder(bits)
    out = []
    while len(out) < limit:
        if stop_bit is not None and c.pos - 32 > stop_bit:
            break
        if not c.code(m.p_more()):
            break
        m.more += 1
        out.append(m.code_byte(c))
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
# 3. Reed-Solomon style erasure code across strands (systematic Cauchy code, GF(256))
# ----------------------------------------------------------------------------------
_EXP, _LOG = [0] * 512, [0] * 256
_v = 1
for _i in range(255):
    _EXP[_i] = _v
    _LOG[_v] = _i
    _v <<= 1
    if _v & 0x100:
        _v ^= 0x11D
for _i in range(255, 512):
    _EXP[_i] = _EXP[_i - 255]


def _gmul(a, b):
    return 0 if a == 0 or b == 0 else _EXP[_LOG[a] + _LOG[b]]


def _ginv(a):
    return _EXP[255 - _LOG[a]]


def _coef(j, i):
    """Coefficient of data strand i in parity strand j: 1 / (x_j + y_i), Cauchy."""
    return _ginv((255 - j) ^ i)


def _parity(data, r):
    K = len(data[0])
    out = []
    for j in range(r):
        p = [0] * K
        for i, d in enumerate(data):
            c = _coef(j, i)
            lc = _LOG[c]
            for t in range(K):
                if d[t]:
                    p[t] ^= _EXP[lc + _LOG[d[t]]]
        out.append(p)
    return out


def _solve_inverse(A):
    """Invert a square matrix over GF(256) (Gauss-Jordan)."""
    m = len(A)
    M = [row[:] + [int(i == j) for j in range(m)] for i, row in enumerate(A)]
    for col in range(m):
        piv = next(r for r in range(col, m) if M[r][col])
        M[col], M[piv] = M[piv], M[col]
        inv = _ginv(M[col][col])
        M[col] = [_gmul(inv, x) for x in M[col]]
        for r in range(m):
            if r != col and M[r][col]:
                f = M[r][col]
                M[r] = [x ^ _gmul(f, y) for x, y in zip(M[r], M[col])]
    return [row[m:] for row in M]


def _recover(data, parity, n):
    """data: {i: bytes-list} for received data strands; parity: {j: list}.
    Returns the full list of n data payloads, or None if not enough survived."""
    missing = [i for i in range(n) if i not in data]
    if not missing:
        return [data[i] for i in range(n)]
    use = sorted(parity)[:len(missing)]
    if len(use) < len(missing):
        return None
    K = len(next(iter(parity.values())))
    rhs = []
    for j in use:
        v = parity[j][:]
        for i, d in data.items():
            c = _coef(j, i)
            for t in range(K):
                if d[t]:
                    v[t] ^= _gmul(c, d[t])
        rhs.append(v)
    Ainv = _solve_inverse([[_coef(j, i) for i in missing] for j in use])
    for a, i in enumerate(missing):
        row = [0] * K
        for b in range(len(use)):
            c = Ainv[a][b]
            if c:
                for t in range(K):
                    if rhs[b][t]:
                        row[t] ^= _gmul(c, rhs[b][t])
        data[i] = row
    return [data[i] for i in range(n)]


# ----------------------------------------------------------------------------------
# 4. Strand framing: [index 8 | n_data 8 | K payload bytes | pad | CRC-8]
# ----------------------------------------------------------------------------------
def _crc(bits, salt=0):
    r = 0x5A ^ salt
    for b in bits:
        r ^= b << 7
        r = ((r << 1) ^ 0x07) & 0xFF if r & 0x80 else (r << 1) & 0xFF
    return [(r >> (CRC_BITS - 1 - i)) & 1 for i in range(CRC_BITS)]


def _to_int(bits):
    x = 0
    for b in bits:
        x = x * 2 + b
    return x


def _to_bits(x, n):
    return [(x >> (n - 1 - i)) & 1 for i in range(n)]


def _K(L):
    """Payload bytes per strand of length L."""
    return (capacity(L) - HDR_BITS - CRC_BITS) // 8


def _frame(index, n, payload, L):
    cap = capacity(L)
    bits = _to_bits(index, 8) + _to_bits(n, 8)
    for byte in payload:
        bits += _to_bits(byte, 8)
    bits += [0] * (cap - CRC_BITS - len(bits))
    bits += _crc(bits, L & 0xFF)
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
    if _crc(bits[:-CRC_BITS], L & 0xFF) != bits[-CRC_BITS:]:
        return None
    index, n = _to_int(bits[:8]), _to_int(bits[8:16])
    K = _K(L)
    payload = [_to_int(bits[16 + 8 * t:24 + 8 * t]) for t in range(K)]
    return index, n, payload


# ----------------------------------------------------------------------------------
# 5. Choosing strand length L and number of parity strands r for each passage
# ----------------------------------------------------------------------------------
def _binom_tail_ok(N, r, q):
    """P(at most r of N strands fail), each failing independently with prob q."""
    p, term = 0.0, (1 - q) ** N
    for k in range(r + 1):
        p += term
        term = term * (N - k) / (k + 1) * q / (1 - q) if q < 1 else 0
    return min(p, 1.0)


def _plan(nbits, nchars):
    best = None
    for L in L_RANGE:
        K = _K(L)
        if K < 4:
            continue
        n = max(1, -(-nbits // (8 * K)))
        if n > 200:
            continue
        for r in range(0, min(12, 255 - n) + 1):
            N = n + r
            # reads per surviving strand, and expected full-length reads among them
            lam = 4 * nchars / (N * (1 - P_LOST)) * (1 - P_BREAK) ** (L - 1)
            # P(fewer than MIN_FULL full reads): too few to find L and vote reliably
            term, p_few = math.exp(-lam), 0.0
            for k in range(MIN_FULL):
                p_few += term
                term *= lam / (k + 1)
            q = P_LOST + (1 - P_LOST) * p_few
            ok = _binom_tail_ok(N, r, q)
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
    strands = [_frame(i, n, d, L) for i, d in enumerate(data)]
    strands += [_frame(n + j, n, p, L) for j, p in enumerate(_parity(data, r))]
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
    # CRC-valid strands (weighted by how many reads support them).
    choice = None
    for L in _candidate_lengths(reads):
        best, nvotes = {}, collections.Counter()
        for cons, support in _consensus_candidates(reads, L):
            got = _unframe(cons)
            if got is None:
                continue
            i, n, payload = got
            nvotes[n] += support
            if i not in best or support > best[i][0]:
                best[i] = (support, payload)
        total = sum(v[0] for v in best.values())
        if best and (choice is None or total > choice[0]):
            choice = (total, L, best, nvotes)
    if choice is None:
        return ""
    _, L, best, nvotes = choice
    n = nvotes.most_common(1)[0][0]
    data = {i: v[1] for i, v in best.items() if i < n}
    parity = {i - n: v[1] for i, v in best.items() if i >= n}
    full = _recover(dict(data), parity, n)
    K = _K(L)
    stop = None
    if full is None:                        # too many losses: decode up to the first gap
        full = []
        for i in range(n):
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
