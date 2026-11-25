import collections
import heapq
import math
import random

# ----- compression: a fixed Huffman code for English -----
# Approximate counts per 100,000 characters of English prose. Built into the code,
# so encode and decode share the same table and nothing has to be sent.
_FREQ = {
    " ": 16500, "e": 9300, "t": 6800, "a": 6000, "o": 5700, "i": 5300, "n": 5300,
    "s": 4900, "r": 4600, "h": 4300, "l": 3100, "d": 3000, "c": 2200, "u": 2100,
    "m": 1900, "f": 1700, "p": 1500, "g": 1500, "w": 1500, "y": 1400, "b": 1100,
    ",": 1000, ".": 900, "v": 750, "k": 500, "\n": 300, "T": 300, "I": 250,
    "A": 200, "'": 200, '"': 200, "-": 150, "S": 150, "x": 130, "H": 100, "W": 100,
    "M": 100, "B": 100, "C": 100, "j": 100, "q": 80, "P": 80, "z": 60, "D": 60,
    "E": 60, "F": 60, "N": 60, "O": 60, "R": 60, "L": 50, "G": 50, "?": 50,
    ";": 30, ":": 30, "!": 30, "J": 30, "Y": 30, "0": 30, "1": 30, "2": 30, "3": 30,
    "4": 30, "5": 30, "6": 30, "7": 30, "8": 30, "9": 30, "K": 20, "U": 20, "V": 20,
    "(": 20, ")": 20, "Q": 5, "X": 5, "Z": 5,
}
for _c in "#$%&*+/<=>@[\\]^_`{|}~\t":
    _FREQ[_c] = 5
_ESC = None             # escape symbol: any character not in the table
_FREQ[_ESC] = 2

def _build_code():
    """Huffman: repeatedly merge the two rarest symbols. Same table -> same code every time."""
    heap = [(f, i, sym) for i, (sym, f) in enumerate(_FREQ.items())]
    heapq.heapify(heap)
    n = len(heap)
    while len(heap) > 1:
        f1, _, a = heapq.heappop(heap)
        f2, _, b = heapq.heappop(heap)
        heapq.heappush(heap, (f1 + f2, n, (a, b)))
        n += 1
    codes = {}
    def walk(node, prefix):
        if isinstance(node, tuple):
            walk(node[0], prefix + "0")
            walk(node[1], prefix + "1")
        else:
            codes[node] = prefix
    root = heap[0][2]
    walk(root, "")
    return codes, root

_CODES, _ROOT = _build_code()

def compress(text):
    bits = [format(len(text), "024b")]                  # 24-bit character count
    for ch in text:
        if ch in _CODES:
            bits.append(_CODES[ch])
        else:                                            # rare character: escape + 21-bit code point
            bits.append(_CODES[_ESC] + format(ord(ch), "021b"))
    bits = "".join(bits)
    bits += "0" * (-len(bits) % 8)
    return int(bits, 2).to_bytes(len(bits) // 8, "big") if bits else b""

def decompress(data):
    bits = "".join(format(b, "08b") for b in data)
    n, pos, out = int(bits[:24], 2), 24, []
    while len(out) < n and pos < len(bits):
        node = _ROOT
        while isinstance(node, tuple) and pos < len(bits):
            node = node[int(bits[pos])]
            pos += 1
        if isinstance(node, tuple):
            break
        if node is _ESC:
            out.append(chr(min(int(bits[pos:pos + 21] or "0", 2), 0x10FFFF)))
            pos += 21
        else:
            out.append(node)
    return "".join(out)

# ----- settings -----
L = 150            # bases per strand
INDEX_BITS = 12    # page number (up to 4096 pages)
COUNT_BITS = 12    # number of data pages, stored on every page
BASES = "ACGT"

# ----- constrained code: never write GG or CC -----
def _allowed(prev):
    """Bases allowed after `prev`: no G after G, no C after C."""
    return [b for b in BASES if not (b == prev and b in "GC")]

_count_cache = {}
def _count(n, prev):
    """How many valid strings of length n can follow the base `prev`."""
    if n == 0:
        return 1
    key = (n, prev)
    if key not in _count_cache:
        _count_cache[key] = sum(_count(n - 1, b) for b in _allowed(prev))
    return _count_cache[key]

def _int_to_strand(x, n):
    """Turn a number into the x-th valid strand of length n (enumerative coding)."""
    out, prev = [], None
    for pos in range(n):
        for b in _allowed(prev):
            c = _count(n - pos - 1, b)
            if x < c:
                out.append(b); prev = b
                break
            x -= c
    return "".join(out)

def _strand_to_int(s):
    """Reverse of _int_to_strand."""
    x, prev, n = 0, None, len(s)
    for pos, ch in enumerate(s):
        for b in _allowed(prev):
            if b == ch:
                break
            x += _count(n - pos - 1, b)
        prev = ch
    return x

STRAND_BITS = _count(L, None).bit_length() - 1      # bits that always fit in one strand (275)
PAYLOAD_BITS = STRAND_BITS - INDEX_BITS - COUNT_BITS # 251
PAYLOAD_BYTES = PAYLOAD_BITS // 8                    # 31

# ----- v3: XOR parity pages that can rebuild lost pages -----
# Parity checks (Lecture 11) applied to whole pages instead of single bits: a parity page is
# the XOR of a group of data pages, so if exactly one page of the group is lost, XOR-ing the
# parity page with the pages that did arrive gives it back.
#
# The groups overlap the way a Hamming code's do. With r parity pages, every data page gets a
# label: an r-bit number with at least two 1s (Hamming skips the powers of two for the same
# reason). Data page i belongs to group j when bit j of its label is 1, so every data page is
# protected by at least two parity pages. Labels with two 1s are used first, then three, ...
#
# Decoding repeats one step until nothing changes: find a group missing exactly one page and
# rebuild that page. Rebuilding one page can leave another group missing only one, and so on.
P_LOST = 0.035      # chance a whole page never comes back (measured on the real channel)

def _labels(r):
    """Labels for the data pages: r-bit numbers with at least two 1s, fewest 1s first."""
    if r == 1:
        return [1]                      # one parity page: the XOR of every data page
    return sorted((x for x in range(1, 2 ** r) if bin(x).count("1") >= 2),
                  key=lambda x: (bin(x).count("1"), x))

def _groups(m, r):
    """groups[j] = the data pages that parity page j covers."""
    if r == 0:
        return []
    labels = _labels(r)
    groups = [[] for _ in range(r)]
    for i in range(m):
        label = labels[i % len(labels)]
        for j in range(r):
            if label >> j & 1:
                groups[j].append(i)
    return groups

def _xor(pages, K):
    out = [0] * K
    for p in pages:
        for t in range(K):
            out[t] ^= p[t]
    return out

def _make_spares(pages, r):
    """pages: list of byte lists. Returns the r parity pages."""
    K = len(pages[0])
    return [_xor([pages[i] for i in g], K) for g in _groups(len(pages), r)]

def _rebuild(have, spares, m, r):
    """have: {i: page} for the data pages that arrived; spares: {j: page} for the parity pages
    that arrived. Rebuilds what it can; pages it cannot rebuild are filled with zeros."""
    K = PAYLOAD_BYTES
    pages = dict(have)
    groups = _groups(m, r)
    progress = True
    while progress and len(pages) < m:
        progress = False
        for j, g in enumerate(groups):
            if j not in spares:
                continue
            missing = [i for i in g if i not in pages]
            if len(missing) == 1:       # this group is missing exactly one page: rebuild it
                pages[missing[0]] = _xor([spares[j]] + [pages[i] for i in g if i != missing[0]], K)
                progress = True
    return [pages.get(i, [0] * K) for i in range(m)]

def _recovers(m, groups, lost):
    """True if the peeling decoder gets back every data page when the pages in `lost` are lost
    (data pages are 0..m-1, parity page j is m+j). groups: list of sets."""
    missing = {i for i in lost if i < m}
    live = [g for j, g in enumerate(groups) if m + j not in lost]
    progress = True
    while missing and progress:
        progress = False
        for g in live:
            gone = g & missing
            if len(gone) == 1:
                missing -= gone
                progress = True
    return not missing

def _p_success(m, r):
    """P(every data page is recovered) when each of the m + r pages is lost with probability
    P_LOST. Law of total probability over k = number of lost pages (binomial): for each k,
    the fraction of k-page losses the decoder survives, counted exactly for k <= 2 and
    estimated by Monte Carlo for k >= 3. The seed is fixed, so encode and decode get the same
    answer."""
    N = m + r
    groups = [set(g) for g in _groups(m, r)]
    rng = random.Random(1000 * m + r)
    def draw(k):                        # k distinct pages, using only rng.random()
        s = set()
        while len(s) < k:
            s.add(int(rng.random() * N))
        return s
    total = 0.0
    for k in range(0, min(N, 8) + 1):
        pk = math.comb(N, k) * P_LOST ** k * (1 - P_LOST) ** (N - k)
        if k == 0:
            frac = 1.0
        elif k == 1:
            frac = sum(_recovers(m, groups, {a}) for a in range(N)) / N
        elif k == 2 and N * (N - 1) // 2 <= 1000:
            frac = sum(_recovers(m, groups, {a, b}) for a in range(N) for b in range(a + 1, N)) / (N * (N - 1) // 2)
        else:
            trials = 200
            frac = sum(_recovers(m, groups, draw(k)) for _ in range(trials)) / trials
        total += pk * frac
    return total

_spares_cache = {}
def _spares_needed(m):
    """How many parity pages to add to m data pages: the number that gives the best
    expected density x chance of success. Tries r = 1, 3, 4, ... and stops once four values
    in a row have not beaten the best so far."""
    if m not in _spares_cache:
        best_r, best = 0, (1 - P_LOST) ** m / m          # r = 0: no parity at all
        worse = 0
        for r in [1] + list(range(3, 60)):
            value = _p_success(m, r) / (m + r)
            if value > best:
                best_r, best, worse = r, value, 0
            elif r >= 3:
                worse += 1
                if worse == 4:
                    break
        _spares_cache[m] = best_r
    return _spares_cache[m]

# ----- encoder and decoder -----
def _layout(n):
    """The data pages and the parity-page numbers. XOR has no size limit (unlike GF(256)),
    so all n data pages form one block."""
    return [(0, n, list(range(n, n + _spares_needed(n))))] if n else []

def encode(text):
    data = compress(text)
    data += bytes(-len(data) % PAYLOAD_BYTES)
    pages = [list(data[p:p + PAYLOAD_BYTES]) for p in range(0, len(data), PAYLOAD_BYTES)]
    n = len(pages)
    for first, m, spare_ids in _layout(n):              # v3: add the parity pages
        pages += _make_spares(pages[first:first + m], len(spare_ids))
    strands = []
    for i, page in enumerate(pages):
        x = (((i << COUNT_BITS) | n) << PAYLOAD_BITS) | int.from_bytes(bytes(page), "big")
        strands.append(_int_to_strand(x, L))
    return strands

def decode(reads):
    votes = collections.defaultdict(collections.Counter)
    n_votes = collections.Counter()
    for r in reads:
        if len(r) == L and "GG" not in r and "CC" not in r:
            x = _strand_to_int(r)
            payload = x & ((1 << PAYLOAD_BITS) - 1)
            n = (x >> PAYLOAD_BITS) & ((1 << COUNT_BITS) - 1)
            i = x >> (PAYLOAD_BITS + COUNT_BITS)
            if payload < (1 << (8 * PAYLOAD_BYTES)) and n > 0:   # skip damaged reads
                votes[i][payload] += 1
                n_votes[n] += 1
    if not votes:
        return ""
    n = n_votes.most_common(1)[0][0]                     # every page carries n, so take a vote
    page = {i: list(votes[i].most_common(1)[0][0].to_bytes(PAYLOAD_BYTES, "big"))
            for i in votes if sum(votes[i].values()) >= 2}
    data = []
    for first, m, spare_ids in _layout(n):
        have = {i - first: page[i] for i in range(first, first + m) if i in page}
        spares = {j: page[s] for j, s in enumerate(spare_ids) if s in page}
        for p in _rebuild(have, spares, m, len(spare_ids)):   # v3: rebuild lost pages if possible
            data += p
    return decompress(bytes(data))
