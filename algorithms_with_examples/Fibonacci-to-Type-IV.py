from functools import reduce
import copy
import random
import itertools

def Monomial(n):
    A = list(range(n))
    M = []
    for i in range(1, 2 ** n):
        Mtemp1 = bin(i)[2:].zfill(n)
        Mtemp2 = list(Mtemp1)
        Mtemp3 = [int(s) for s in Mtemp2]
        Mtemp4 = []
        for j in range(n):
            if Mtemp3[j] == 1:
                Mtemp4.append(A[j])
        M.append(Mtemp4)
    M.remove([0])
    return M


def Fib_NLFSR(R, FF, ZF, N0):
    N = copy.deepcopy(N0)
    States = []
    output = []
    for i in range(R):
        if i <= 2**n:
            States.append(N)
            #print('NFib=',N)
        # First compute the outputs
        Z = copy.deepcopy(ZF)
        for j in range(len(Z)):
            Z[j] = list(map(lambda x: N[x], Z[j]))
        for j in range(len(Z)):
            Z[j] = reduce(lambda p, q: p & q, Z[j])
        #for j in range(len(Z)):
            #Z = reduce(lambda p, q: p ^ q, Z)
        Z = reduce(lambda p, q: p ^ q, Z)
        output.append(Z)

        # Second compute the internal states for next round
        F = copy.deepcopy(FF)
        for j in range(len(F)):
            for k in range(len(F[j])):
                F[j][k] = list(map(lambda x: N[x], F[j][k]))
        for j in range(len(F)):
            for k in range(len(F[j])):
                F[j][k] = reduce(lambda p, q: p & q, F[j][k])
        for j in range(len(F)):
            F[j] = reduce(lambda p, q: p ^ q, F[j])
        N = copy.deepcopy(F)
    return output, States


def Gal_NLFSR(R, FF, ZF, N0):
    N = copy.deepcopy(N0)
    States = []
    output = []
    for i in range(R):
        if i <= 2**n:
            States.append(N)
            #print('N=',N)
        # First compute the outputs
        Z = copy.deepcopy(ZF)
        for j in range(len(Z)):
            Z[j] = list(map(lambda x: N[x], Z[j]))
        for j in range(len(Z)):
            Z[j] = reduce(lambda p, q: p & q, Z[j])
        #for j in range(len(Z)):
            #Z = reduce(lambda p, q: p ^ q, Z)
        Z = reduce(lambda p, q: p ^ q, Z)
        output.append(Z)

        # Second compute the internal states for next round
        F = copy.deepcopy(FF)
        for j in range(len(F)):
            for k in range(len(F[j])):
                F[j][k] = list(map(lambda x: N[x], F[j][k]))
        for j in range(len(F)):
            for k in range(len(F[j])):
                F[j][k] = reduce(lambda p, q: p & q, F[j][k])
        for j in range(len(F)):
            F[j] = reduce(lambda p, q: p ^ q, F[j])
        N = copy.deepcopy(F)
    return output, States


def FuncSort(Fs): # Remove all duplicate product terms in each feedback function in set Fs. Example, Fs = [[[1]], [[2]], [[3]], [[4]], [[0], [2], [3], [1, 2], [2, 3], [1, 2]]]
    Fs_Sort = []
    for i in range(len(Fs)): # For each feedback function in Fs
        tempx = Fs[i]
        tempy = []
        for sublist in tempx:
            if sublist not in tempy:
                tempy.append(sublist)
            else:
                tempy.remove(sublist)
        Fs_Sort.append(tempy)
    return Fs_Sort

def Compu_B_Backw(MBackw, n):
    BBackw = []
    for i in range(len(MBackw)):
        m = copy.deepcopy(MBackw[i])
        b = n - 1 - max(m) - 1
        BBackw.append(b)
    return BBackw

def Compen_List_Backw(MBackw, BBackw, n):
    CBackw = [-1] * n
    for i in range(len(MBackw)):
        C = [-1] * n
        m = copy.deepcopy(MBackw[i])
        b = copy.deepcopy(BBackw[i])
        for j in range(b+1):
            CTemp = list(map(lambda  x: x + j, m))
            C[j] = copy.deepcopy(CTemp)
        for j in range(n):
            if CBackw[j] == -1 and C[j] != -1:
                CBackw[j] = copy.deepcopy([C[j]])
            elif CBackw[j] != -1 and C[j] != -1:
                CBackw[j].append(C[j])
    return CBackw

def Shifting_Backw(FF, MBackw, BBackw, n):
    FFShitedBackw = copy.deepcopy(FF)
    for i in range(len(MBackw)):
        MTemp = list(map(lambda x: x + BBackw[i] + 1, MBackw[i]))
        FFShitedBackw[n-1].append(MBackw[i])
        FFShitedBackw[BBackw[i]].append(MTemp)
    FFShitedBackw = FuncSort(FFShitedBackw)
    return FFShitedBackw

def Compensation_Backw(FF, Z, CBackw, n):
    FFCompenBackw = copy.deepcopy(FF)
    FFCompenBackw.append(Z)
    for i in range(n):
        if CBackw[i] != -1:
            for j in range(n+1):
                FFj = copy.deepcopy(FFCompenBackw[j])
                for k in range(len(FFCompenBackw[j])):
                    if j < n and FFCompenBackw[j][k] == [(j + 1) % n]:
                        continue
                    FFjk = copy.deepcopy(FFCompenBackw[j][k]) # FFjk = [2, 3]
                    if i in FFjk:
                        if FFjk in FFj:
                            FFj.remove(FFjk)
                        if len(FFCompenBackw[j][k]) == 1:
                            temp = copy.deepcopy(CBackw[i]+[[i]])
                        else:
                            FFjk = copy.deepcopy(FFCompenBackw[j][k])
                            FFjk.remove(i)
                            temp = list(map(lambda x: list(set(x+FFjk)), CBackw[i]+[[i]]))
                        FFj = FFj + temp
                FFCompenBackw[j] = copy.deepcopy(FFj)
            FFCompenBackw = FuncSort(FFCompenBackw)
            ZCompenBackw = copy.deepcopy(FFCompenBackw[-1])
    return FFCompenBackw[:-1], ZCompenBackw

def Compute_IV(N0, C, n):
    N0Gal = copy.deepcopy(N0)
    CTemp = copy.deepcopy(C)
    for i in range(n):
        if CTemp[i] != -1:
            gi = CTemp[i] + [[i]]
            #print('gi=', gi)
            for j in range(len(gi)):
                gi[j] = list(map(lambda x: N0[x], gi[j]))
            #print('gi1=',gi)
            for j in range(len(gi)):
                gi[j] = reduce(lambda p, q: p & q, gi[j])
            #print('gi2=', gi)
            gi = reduce(lambda p, q: p ^ q, gi)
            #print('gi3=', gi)
            N0Gal[i] = gi
    return N0Gal

n = 7
n = 4
R = 16
N0 = []  # Generate the initial value for the NLFSR
for i in range(n):
    N0.append(random.randint(0, 1))
#N0 = [1, 1, 1, 0]
print('N0=', N0)

"""
M = Monomial(n)
MTemp = copy.deepcopy(M)
for i in range(len(MTemp)):
    if MTemp[i][0] == 0 or len(MTemp[i]) > 2:
        M.remove(MTemp[i])
print('M=', M)
length = min(n, 3)
F = random.sample(M, length)
F.insert(0, [0])
"""
F = [[0], [3], [3, 5], [1, 2]]
F = [[0], [1]]
#F = [[0], [1], [2], [1, 2]]
FF = []
for j in range(n - 1):
    f = [[j + 1]]
    FF.append(f)
FF.append(F)
print('FF=', FF)

Z = [[0]]

M = copy.deepcopy(F[1:])
#M = [[1], [1, 2]]
B = Compu_B_Backw(M, n)
#B = [0]
print('B =', B)

BAll = []
for i in range(len(B)):
    BTemp = list(range(0, B[i]+1, 1)) #n-1
    BAll.append(BTemp)
BAll = list(itertools.product(*BAll))
BAll = [list(e) for e in BAll]
print('BAll=', BAll)
print('lenBAll=', len(BAll))

for i in range(len(BAll)):
    FFFib = copy.deepcopy(FF)
    ZFib = copy.deepcopy(Z)
    N0Fib = copy.deepcopy(N0)
    MBackw = copy.deepcopy(M)
    BBackw = copy.deepcopy(BAll[i])
    print('FFFib =', FFFib)
    print('ZFib =', ZFib)
    print('N0Fib =', N0Fib)

    CBackw = Compen_List_Backw(MBackw, BBackw, n)
    print('CBackw =', CBackw)

    FFShiftedBackw = Shifting_Backw(FFFib, MBackw, BBackw, n)
    FFGal, ZGal = Compensation_Backw(FFShiftedBackw, ZFib, CBackw, n)
    N0Gal = Compute_IV(N0, CBackw, n)
    print('FFGal =', FFGal)
    print('ZGal =', ZGal)
    print('N0Gal =', N0Gal)

    OutFib, StateFib = Fib_NLFSR(R, FFFib, ZFib, N0Fib)
    OutGal, StateGal = Gal_NLFSR(R, FFGal, ZGal, N0Gal)
    ZGal2 = [[0]]
    OutGal2, StateGal2 = Gal_NLFSR(R, FFGal, ZGal2, N0Gal)
    if OutFib == OutGal:
        print('Success!')
        print('StateFib =', StateFib)
        print('OutFib =', OutFib)
        print('StateGal =', StateGal)
        print('OutGal2 =', OutGal2)
    else:
        print('Fail!!!')
        break







