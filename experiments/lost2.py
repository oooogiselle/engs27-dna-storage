import requests, random, collections, pickle
U="https://dna.cedar.host.dartmouth.edu"
def reads(strands,n):
    got=[]
    while len(got)<n:
        k=min(n-len(got),5000)
        got+=requests.post(U+"/api/reads",json={"key":"","strands":strands,"n":k},timeout=120).json()["reads"]
    return got
def count(S,R,k=12):
    idx=collections.defaultdict(set)
    for j,x in enumerate(S):
        for i in range(len(x)-k+1): idx[x[i:i+k]].add(j)
    c=collections.Counter()
    for r in R:
        if len(r)>=k:
            h=idx.get(r[:k],set())
            if len(h)==1: c[next(iter(h))]+=1
    return c
S,c=pickle.load(open("lost.pkl","rb"))
lost=[S[j] for j in range(len(S)) if c[j]==0][:60]
kept=[S[j] for j in range(len(S)) if c[j]>=15][:240]
# 1. determinism: resend 60 previously-lost with 240 previously-kept, twice
for rep in range(2):
    T=lost+kept; R=reads(T,15000); cc=count(T,R)
    print("rep",rep,"prev-lost still lost:",sum(cc[j]==0 for j in range(60)),"/60  prev-kept now lost:",sum(cc[j]==0 for j in range(60,300)),"/240",flush=True)
# 2. length dependence
rng=random.Random(7)
def gen(L):
    s=[rng.choice("ACGT")]
    while len(s)<L:
        b=rng.choice("ACGT")
        if b==s[-1] and b in "GC": continue
        s.append(b)
    return "".join(s)
for L in (40,100,200,400):
    n=min(800,90000//L); T=[gen(L) for _ in range(n)]
    R=reads(T,20*n); cc=count(T,R)
    print("L",L,"n",n,"lost frac",round(sum(cc[j]==0 for j in range(n))/n,3),"median reads",sorted(cc[j] for j in range(n))[n//2],flush=True)
