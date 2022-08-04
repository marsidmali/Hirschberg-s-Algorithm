import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-t", action="store_true", help="print i,j")
parser.add_argument("gap", type=int, help="gap penalty")
parser.add_argument("match", type=int, help="match")
parser.add_argument("diff", type=int, help="diff")
parser.add_argument("seq1",type=str, help="1st sequence to be aligned")
parser.add_argument("seq2",type=str, help="2nd sequence to be aligned")

args = parser.parse_args()
gap = args.gap
match = args.match
diff = args.diff
t = args.t
A = args.seq1
B = args.seq2


def Compare(A, B): 
    if A == B:
        return match
    else:
        return diff


def NeedlemanWunsch(A, B): 
    F = [[0 for x in range(len(B) + 1)] for x in range(len(A) + 1)]  
    for i in range(len(A)+1):  
        F[i][0] = i * gap
    for j in range(len(B)+1):  
        F[0][j] = j * gap
    for i in range(1, len(A) + 1):  
        for j in range(1, len(B) + 1):
            md = Compare(A[i - 1], B[j - 1])
            F[i][j] = max(F[i - 1][j] + gap, F[i][j - 1] + gap, F[i - 1][j - 1] + md)
    global WW, ZZ
    WW = []
    ZZ = []
    WW, ZZ = EnumerateAlignments(A, B, F, "", "")
    return WW, ZZ


def EnumerateAlignments(A, B, F, W, Z): 
    i = len(A)
    j = len(B)

    if i == 0 and j == 0:
        WW.append(W)
        ZZ.append(Z)
    if i > 0 and j > 0:
        md = Compare(A[i - 1], B[j - 1])
        if F[i][j] == F[i - 1][j - 1] + md:  
            EnumerateAlignments(A[:i - 1], B[:j - 1], F, A[i - 1] + W, B[j - 1] + Z)
    if i > 0 and F[i][j] == F[i - 1][j] + gap:  
        EnumerateAlignments(A[:i - 1], B, F, A[i - 1] + W, "-" + Z)
    if j > 0 and F[i][j] == F[i][j - 1] + gap: 
        EnumerateAlignments(A, B[:j - 1], F, "-" + W, B[j - 1] + Z)
    return WW, ZZ


def ComputeAlignmentScore(A, B):
    L = []
    K = []
    for i in range(len(B)+1):
        L.append(i * gap)
        K.append(0)
    for i in range(1, len(A) + 1):
        L, K = K, L  
        L[0] = i * gap  
        for j in range(1, len(B) + 1):  
            md = Compare(A[i - 1], B[j - 1])
            L[j] = max(L[j - 1] + gap, K[j] + gap, K[j - 1] + md)  
    return L


def UpdateAlignments(WW, ZZ, WWl,WWr,ZZl,ZZr): 
    for i, j in zip(WWl, ZZl):  
        for k, l in zip(WWr, ZZr):
            WW.append(i + k)
            ZZ.append(j + l)


def Hirschberg(A,B): 
    WW,ZZ = [],[]
    if len(A) == 0:  
        WW.append('-' * len(B))
        ZZ.append(B)
    elif len(B) == 0:  
        WW.append(A)
        ZZ.append('-' * len(A))
    elif len(A) == 1 or len(B) == 1:
        WW,ZZ = NeedlemanWunsch(A,B)
    else:
        i = len(A) // 2
        Sl = ComputeAlignmentScore(A[:i], B)   
        Sr = ComputeAlignmentScore(A[i:][::-1],B[::-1])  
        S = [m + n for m,n in zip(Sl, Sr[::-1])] 
        J = [x for x in range(len(S)) if S[x] == max(S)] 
        WW = []
        ZZ = []
        for j in J:  
            if t:
                print(str(i)+",",j)
            WWl, ZZl = Hirschberg(A[:i], B[:j])
            WWr, ZZr = Hirschberg(A[i:len(A)], B[j:len(B)])
            UpdateAlignments(WW, ZZ, WWl, WWr, ZZl, ZZr)
    return WW, ZZ


def print_alignments()
    WW, ZZ = Hirschberg(A, B)
    a = []
    for i in zip(WW, ZZ):
        if i not in a:
            a.append(i)
    for w, z in a:
        print(w)
        print(z +"\n")


print_alignments()
