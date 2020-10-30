"""
The implementation of the Type-I-to-Fibonacci transformation algorithm.
"""

from functools import reduce
#from itertools import permutations
from itertools import product
import copy
#import random
#import itertools

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
    #M.remove([0])
    return M

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
        if Fs[i] == 0:
            break
        tempx = Fs[i]
        tempy = []
        for sublist in tempx:
            if sublist not in tempy:
                tempy.append(sublist)
            else:
                tempy.remove(sublist)
        Fs_Sort.append(tempy)
    return Fs_Sort

def Shift_Single_Gi(FF, t, n):
    FFGal = copy.deepcopy(FF)
    FFGal[t] = [[t + 1]]
    Ci = [0]*n
    for i in range(t + 1):
        Ci[i] = [[(n - (t - i) - 1) % n]]
    for i in range(t - 1, -1, -1):
        CiTemp = copy.deepcopy(Ci[i])
        for j in range(len(CiTemp)):
            if len(FFGal[CiTemp[j][0]]) > 1:
                for k in range(i + 1):
                    Ci[k].append([n - (i - k) - 1])
            #print('Ci =', Ci)
        #Ci = FuncSort(Ci)

    for i in range(t, -1, -1):
        for j in range(1, len(FFGal[-1])):
            Temp = []
            if i in FFGal[-1][j]:
                temp = copy.deepcopy(FFGal[-1][j])
                temp.remove(i)
                Temp.append([temp])
                Temp = Temp + [Ci[i]]
                #print('Temp =', Temp)
                Temp = list(product(*Temp)) #combine elements of lists in Temp3
                Temp = [list(elem) for elem in Temp]
                for k in range(len(Temp)): #convert list of lists into a list
                    Temp[k] = [item for sublist in Temp[k] for item in sublist]
                    Temp[k] = list(set(Temp[k]))
                #print('Temp1 =', Temp)
                #FFGal[-1].remove(FFGal[-1][j])
                FFGal[-1] += Temp
                #print('FFGal =', FFGal)
    FFGal[-1] += Ci[0]
    for i in range(n - 1):
        if len(FFGal[i]) > 1:
            FFGal[i] = [[i + 1]] + FFGal[-1]
    #FFGal = FuncSort(FFGal)
    return Ci, FFGal


n = 5
R = 1000

#FF1 = [[[1]], [[2]], [[3], [0], [1, 2]], [[4], [0], [1, 2]], [[5], [0], [1, 2]], [[0], [1, 2]]]
#FF1 = [[[1]], [[2], [0], [1, 2]], [[3], [0], [1, 2]], [[4], [0], [1, 2]], [[0], [1, 2]]]
FF1 = [[[1], [0], [1, 2]], [[2], [0], [1, 2]], [[3]], [[4], [0], [1, 2]], [[0], [1, 2]]]
ZF1 = [[4]]
N01 = [0, 1, 1, 0, 1]

C = []
FF2 = copy.deepcopy(FF1)
for i in range(n-1):
    if len(FF2[i]) > 1:
        Ci, FF2 = Shift_Single_Gi(FF2, i, n)
        FF2 = FuncSort(FF2)
        print('Ci =', Ci)
        if C == []:
            C = copy.deepcopy(Ci)
        else:
            for j in range(n):
                if Ci[j] != 0 and C[j] != 0:
                    C[j] = C[j] + Ci[j]
                elif Ci[j] != 0 and C[j] == 0:
                    C[j] = copy.deepcopy(Ci[j])
print('FF2 =', FF2)
print('C =', C)
ZF2 = [[4]]
N02 = [0, 0, 0, 1, 1]

Output1, States1 = Gal_NLFSR(R, FF1, ZF1, N01)
Output2, States2 = Gal_NLFSR(R, FF2, ZF2, N02)
print('Output1 =', Output1)
print('Output2 =', Output2)
if Output1 == Output2:
    print('Equivalent!')
