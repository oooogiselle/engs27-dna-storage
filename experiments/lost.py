import requests, random, collections, pickle, sys
U="https://dna.cedar.host.dartmouth.edu"
def reads(strands,n):
    got=[]
    while len(got)<n:
        k=min(n-len(got),5000)
        r=requests.post(U+"/api/reads",json={"key":"","strands":strands,"n":k},timeout=600)
        got+=r.json()["reads"]
    return got
rng=random.Random(11)
def gen(L,ban):
    s=[rng.choice("ACGT")]
    while len(s)<L:
        c=rng.choice("ACGT")
        if ban and c==s[-1] and c in "GC": continue
        s.append(c)
    return "".join(s)
allS,cnt=[],collections.Counter()
for batch in range(4):
    S=[gen(100, batch%2==1) for _ in range(900)]
    R=reads(S,20000)
    idx=collections.defaultdict(list)
    for j,x in enumerate(S):
        for i in range(0,100-11): idx[x[i:i+12]].append(j)
    c=collections.Counter()
    for r in R:
        if len(r)>=12:
            h=set(idx.get(r[:12],[]))
            if len(h)==1: c[h.pop()]+=1
    base=len(allS)
    for j in range(len(S)): cnt[base+j]=c[j]
    allS+=S
    v=sorted(c[j] for j in range(len(S)))
    print("batch",batch,"ban",batch%2==1,"lowest",v[:15],"median",v[450],flush=True)
pickle.dump((allS,cnt),open("lost.pkl","wb"))
