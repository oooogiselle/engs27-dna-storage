import pickle, collections, random
S,R,loc=pickle.load(open("data.pkl","rb"))
L=200
exact=set(r for r in R if any(r in s for s in S)) if False else None
sidx=collections.defaultdict(list)
K=7
for j,s in enumerate(S):
    for i in range(L-K+1): sidx[s[i:i+K]].append((j,i))
def semi(r,s):
    # semi-global: r fully aligned, s free ends; returns (cost, ops string counts)
    n,m=len(r),len(s)
    D=[[0]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): D[i][0]=i
    for i in range(1,n+1):
        for j in range(1,m+1):
            D[i][j]=min(D[i-1][j]+1,D[i][j-1]+1,D[i-1][j-1]+(r[i-1]!=s[j-1]))
    # free end in s
    jm=min(range(m+1),key=lambda j:D[n][j]); cost=D[n][jm]
    # traceback
    i,j=n,jm; ops=collections.Counter()
    while i>0:
        if j>0 and D[i][j]==D[i-1][j-1]+(r[i-1]!=s[j-1]):
            ops["sub" if r[i-1]!=s[j-1] else "match"]+=1; i-=1; j-=1
        elif D[i][j]==D[i-1][j]+1: ops["ins"]+=1; i-=1   # extra base in read
        else: ops["del"]+=1; j-=1
    return cost,ops
random.seed(0)
exactset={}
bad=[]
for r in R:
    if not any(r in s for s in S[:0]): pass
# determine exact quickly via loc counts is awkward; recompute
allsub=set()
for s in S:
    pass
bad=[r for r in random.sample(R,6000) if not any(r in s for s in S)][:600]
tot=collections.Counter(); costs=collections.Counter(); basestot=0; lens=[]
for r in bad:
    votes=collections.Counter()
    for i in range(0,len(r)-K+1):
        for j,p in sidx.get(r[i:i+K],[]): votes[(j,p-i)]+=1
    if not votes: tot["unplaced"]+=1; continue
    (j,off),_=votes.most_common(1)[0]
    s=S[j][max(0,off-4):off+len(r)+6]
    c,ops=semi(r,s)
    costs[c]+=1; tot.update(ops); basestot+=len(r); lens.append(len(r))
print("bad analysed",len(bad),"unplaced",tot["unplaced"])
print("cost dist",sorted(costs.items()))
print({k:v for k,v in tot.items()}, "per base sub",tot["sub"]/basestot,"ins",tot["ins"]/basestot,"del",tot["del"]/basestot)
# error rate over ALL reads: 
print("frac reads with error (est)",1-sum(1 for r in R[:5000] if any(r in s for s in S))/5000)
