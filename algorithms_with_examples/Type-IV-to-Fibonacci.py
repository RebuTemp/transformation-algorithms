from functools import reduce
import copy
import random
import itertools

def MonomialRange(a, b):    #Generate all the possible combinations of {x_a, ..., x_b-1}
    A = list(range(a, b))
    M = []
    for i in range(1, 2 ** (len(A))):
        Mtemp1 = bin(i)[2:].zfill(len(A))
        Mtemp1 = list(Mtemp1)
        Mtemp1 = [int(e) for e in Mtemp1]
        #print('Mtemp1 = ', Mtemp1)
        Mtemp2 = []
        for j in range(len(A)):
            if Mtemp1[j] == 1:
                Mtemp2.append(A[j])
        M.append(Mtemp2)
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


#Given a single monomial such as x2x3, compute the Compensation set C
def CompenSet(b, M, C): # b is the index of the function; M is a single monomial such as [2, 3]; C is the compensation set for M
    for i in range(b + 1):
        Temp = list(map(lambda x: x - (b - i) - 1, M))
        if C[i] == -1:
            C[i] = []
        C[i].append(Temp)
    return C


#Compensate the feedback functions by the compensation set from right to left
def Compensation(Compen, FFs):
    FFAll = copy.deepcopy(FFs)
    for i in range(len(FFAll)):
        if len(FFAll[i]) > 1:
            CompenStart = i
            break

    for i in range(len(Compen)-1, -1, -1):
        if Compen[i] != -1:
            for j in range(CompenStart, n):
                FFsj = copy.deepcopy(FFAll[j])
                for k in range(len(FFAll[j])):
                    if FFAll[j][k] == [(j + 1) % n]:
                        continue
                    FFsjk = copy.deepcopy(FFAll[j][k]) # FFsjk = [2, 3]
                    if i in FFsjk:
                        if FFsjk in FFsj:
                            FFsj.remove(FFsjk)
                        if len(FFAll[j][k]) == 1:
                            temp = copy.deepcopy(Compen[i])
                        else:
                            FFsjk = copy.deepcopy(FFAll[j][k])
                            FFsjk.remove(i)
                            temp = list(map(lambda x: list(set(x+FFsjk)), Compen[i]))
                        FFsj = FFsj + temp
                FFAll[j] = copy.deepcopy(FFsj)
            FFAll = FuncSort(FFAll)
            #print('FFAll=', FFAll)
    return FFAll


#The entire transformation of the feedback function and the output function
def Transformation(n, FFGal):
    FFTemp = copy.deepcopy(FFGal)
    B = [] #The starting position of the shifted monomials
    MGal_All = [] #The monomials to be shifted
    MFib_All = [] #The monomials after shifting
    C_All = [-1] * n
    for i in range(n - 2, -1, -1): #Consider all the feedback functions except the last bit
        if len(FFTemp[i]) > 1: #After each loop, the feedback functions may changed, so we use FFTemp here!
            B.append(i)
            MGal = copy.deepcopy(FFTemp[i][1:])
            MFib = []
            Ci = [-1] * n
            for j in range(len(MGal)):
                Temp1 = list(map(lambda x: x - i - 1, MGal[j]))
                MFib.append(Temp1) #Compute the indexes of monomials in gi after shifting
                Ci = CompenSet(i, MGal[j], Ci)
            MGal_All.append(MGal)
            MFib_All.append(MFib)
            for j in range(n): #Combine all the compensation sets
                if C_All[j] != -1 and Ci[j] != -1:
                    C_All[j] = C_All[j] + Ci[j]
                elif C_All[j] == -1:
                    C_All[j] = copy.deepcopy(Ci[j])
            FFTemp[i] = [[(i + 1) % n]] #Remove the shifted monomials gi in fi
            for j in range(n):
                if Ci[j] != -1:
                    Ci[j].insert(0, [j])
            FFTemp = Compensation(Ci, FFTemp)
        print('FFTemp = ', FFTemp)
    #print('B = ', B)
    #print('MGal_All = ', MGal_All)
    #print('MFib_All = ', MFib_All)
    for i in range(n):
        if C_All[i] != -1:
            C_All[i].insert(0, [i])
    #print('C_All = ', C_All)
    FFFib = copy.deepcopy(FFTemp)
    for i in range(len(MFib_All)): #Add the shifted monomials to f_n-1
        FFFib[-1] = FFFib[-1] + MFib_All[i]
    #print('FFFib = ', FFFib)
    return FFFib, C_All, B, MGal_All, MFib_All


#Compensate the initial value
def CompenIV(n, Compen, N0Origi):
    N0Trans = copy.deepcopy(N0Origi)
    for i in range(n-1, -1, -1):
        if Compen[i] != -1:
            Temp = copy.deepcopy(Compen[i])
            for j in range(len(Compen[i])):
                Temp[j] = list(map(lambda x: N0Trans[x], Temp[j]))
            for j in range(len(Compen[i])):
                Temp[j] = reduce(lambda p, q: p & q, Temp[j])
            Temp = reduce(lambda p, q: p ^ q, Temp)
            N0Trans[i] = Temp
    return N0Trans




n = 7 #length of the Gal_NLFSR
R = 2000 #The number of clocks we run the Gal_NLFSR
N0Gal = []  # Generate a random initial value for the Gal_NLFSR
for i in range(n):
    N0Gal.append(random.randint(0, 1))
#N0Gal = [1, 1, 0, 0]
print('N0Gal = ', N0Gal)


#Random feedback functions and output function for the Gal_NLFSR
FFGal = []
for i in range(4):
    M = MonomialRange(i + 2, n)
    #print('M = ', M)
    MTemp = copy.deepcopy(M)
    for j in range(len(MTemp)):
        if MTemp[j][0] == 0 or len(MTemp[j]) > 2:
            M.remove(MTemp[j])
    #print('M2 = ', M)
    length = random.randint(0, 2) #Choose the random number of monomials we add to the feedback functions
    length = min(length, len(M))
    f = random.sample(M, length)
    f.insert(0, [(i + 1) % n])
    FFGal.append(f)
for i in range(4, n):
    FFGal.append([[(i + 1) % n]])
FFGal[-1].append([4, 5])
#FFGal =  [[[1], [6, 5]], [[2], [4, 5]], [[3], [4]], [[4]], [[5]], [[6]], [[0], [4, 5]]]
#FFGal = [[[1]], [[2], [4]], [[3]], [[4], [5]], [[5]], [[6]], [[0], [4, 5]]]
#FFGal = [[[1]], [[2], [4], [6]], [[3]], [[4], [5]], [[5]], [[6]], [[0], [4, 5]]]
#FFGal = [[[1]], [[2], [3]], [[3]], [[4], [5]], [[5]], [[6]], [[0], [4, 5]]]
#FFGal =  [[[1], [2, 4]], [[2], [4, 5]], [[3], [4]], [[4]], [[5]], [[6]], [[0], [4, 5]]]
#FFGal = [[[1]], [[2], [4]], [[3]], [[4], [5]], [[5]], [[6]], [[0], [4, 5]]]
FFGal = [[[1], [5, 6], [5]], [[2], [6]], [[3]], [[4]], [[5]], [[6]], [[0], [4, 5]]]
#FFGal = [[[1], [5, 6]], [[2], [6]], [[3]], [[4]], [[5]], [[6]], [[0], [4, 5]]]
print('FFGal = ', FFGal)

z = 0 #Choose the bit to generate the output
ZGal = [[z]]
print('ZGal = ', ZGal)
OutGal, StateGal = Gal_NLFSR(R, FFGal, ZGal, N0Gal) #Run the Gal_NLFSR
print('OutGal = ', OutGal)

#Compensate the feedback functions and the initial value
FFFib, C_All, B, MGal_All, MFib_All = Transformation(n, FFGal)
#FFFib = [[[1]], [[2]], [[3]], [[4]], [[5]], [[6]], [[0], [4, 5], [1]]]
N0Fib = CompenIV(n, C_All, N0Gal)
print('N0Fib = ', N0Fib)
print('FFFib = ', FFFib)
print('MGal_All = ', MGal_All)
print('C_All = ', C_All)

#Construct output function and run the Fib_NLFSR
ZFib_All = copy.deepcopy(C_All)
for i in range(n):
    if ZFib_All[i] == -1:
        ZFib_All[i] = [[i]]
ZFib = ZFib_All[z]
print('ZFib = ', ZFib)

#Run the Fib_NLFSR and compare the output sequences
OutFib, StateFib = Fib_NLFSR(R, FFFib, ZFib, N0Fib)
print('OutFib = ', OutFib)
if OutFib == OutGal:
    print('Yes!')
    #print('MGal_All = ', MGal_All)
    #print('B = ', B)
    #print('MFib_All = ', MFib_All)
    #print('C_All = ', C_All)
    #print('FFGal = ', FFGal)
    #print('ZGal = ', ZGal)
    #print('N0Gal = ', N0Gal)
    #print('OutGal = ', OutGal)
    #print('FFFib = ', FFFib)
    #print('ZFib = ', ZFib)
    #print('N0Fib = ', N0Fib)
    #print('OutFib = ', OutFib)

