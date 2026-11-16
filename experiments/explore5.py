import requests, random, collections, pickle
U="https://dna.cedar.host.dartmouth.edu"
def reads(strands,n):
    got=[]
    while len(got)<n:
        k=min(n-len(got),5000)
        got+=requests.post(U+"/api/reads",json={"key":"","strands":strands,"n":k},timeout=300).json()["reads"]
    return got
rng=random.Random(5)
def rs(L,ban=True):
    s=[rng.choice("ACGT")]
    while len(s)<L:
        c=rng.choice("ACGT")
        if ban and c==s[-1] and c in "GC": continue
        s.append(c)
    return "".join(s)
S0,R0,loc0=pickle.load(open("data.pkl","rb"))
# error rate vs read length on random strands
byl=collections.defaultdict(lambda:[0,0])
for r in R0[:60000]:
    ok=any(r in s for s in S0)
    byl[min(len(r)//10,10)][0]+=1; byl[min(len(r)//10,10)][1]+= (not ok)
print("err frac by len/10:",{k:round(v[1]/v[0],3) for k,v in sorted(byl.items())})
for L in (30,45,60,80,120):
    for ban in (False,True):
        S=[rs(L,ban) for _ in range(60)]
        R=reads(S,10000)
        ln=[len(r) for r in R]
        full=sum(1 for r in R if len(r)==L)/len(R)
        ex=sum(1 for r in R if any(r in s for s in S))/len(R)
        print(f"L={L} ban={ban} mean={sum(ln)/len(ln):.1f} full={full:.3f} exact={ex:.3f}")
# homopolymer / repeated-kmer tests
H={"AAAA":"A"*100,"TTTT":"T"*100,"GGGG":"G"*100,"ACAC":"AC"*50,"AT":"AT"*50}
for k,v in H.items():
    R=reads([v],2000); print(k,"mean len",sum(map(len,R))/len(R))
