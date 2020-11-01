"""
The implementation of the Type-IV-to-Fibonacci transformation algorithm.
"""

import copy
import random
import itertools
from functools import reduce


def Fib_NLFSR(R, FF, ZF, N0):
    """Run Fibonacci NLFSR
    Parameters
    ----------
    R : int
        the number of output bits
    FF : list
        feedback functions of the Fibonacci NLFSR
    ZF : list
        output function of the Fibonacci NLFSR
    N0 : list
        initial state of the Fibonacci NLFSR
    Return
    ----------
    outputs : list
        output bits
    states: list
        internal state for each round
    """

    N = copy.deepcopy(N0)
    states, outputs = [], []
    for i in range(R):
        if i <= 2**n:
            states.append(N)

        # First, compute the outputs
        Z = copy.deepcopy(ZF)
        for j in range(len(Z)):
            Z[j] = list(map(lambda x: N[x], Z[j]))

        for j in range(len(Z)):
            Z[j] = reduce(lambda p, q: p & q, Z[j])
        Z = reduce(lambda p, q: p ^ q, Z)
        outputs.append(Z)

        # Second, compute the internal states for next round
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

    return outputs, states


def Gal_NLFSR(R, FF, ZF, N0):
    """Run Galois NLFSR
    Parameters
    ----------
    R : int
        the number of output bits
    FF : list
        feedback functions of the Galois NLFSR
    ZF : list
        output function of the Galois NLFSR
    N0 : list
        initial state of the Galois NLFSR
    Return
    ----------
    outputs : list
        output bits
    states: list
        internal state for each round
    """

    N = copy.deepcopy(N0)
    states, outputs = [], []
    for i in range(R):
        if i <= 2**n:
            states.append(N)

        # First, compute the outputs
        Z = copy.deepcopy(ZF)
        for j in range(len(Z)):
            Z[j] = list(map(lambda x: N[x], Z[j]))

        for j in range(len(Z)):
            Z[j] = reduce(lambda p, q: p & q, Z[j])
        Z = reduce(lambda p, q: p ^ q, Z)
        outputs.append(Z)

        # Second, compute the internal states for next round
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

    return outputs, states


def sort_function(FF):
    """Remove all duplicate monomials in the feedback functions FF.
    """

    FF_Sort = []
    for i in range(len(FF)):  # For each feedback function in FF
        if type(FF[i]) == list:
            tempx = FF[i]
            tempy = []
            for sublist in tempx:
                if sublist not in tempy:
                    tempy.append(sublist)
                else:
                    tempy.remove(sublist)
            FF_Sort.append(tempy)
        else:
            FF_Sort.append(FF[i])

    return FF_Sort


def compensation_list(b, M, C):
    """Construct compensation list
    Parameters
    ----------
    b : int
        the index of the feedback function
    M : list
        a monomial
    C : list
        an initial list for M
    Return
    ----------
    C : list
        compensation list
    """

    for i in range(b + 1):
        Temp = list(map(lambda x: x - (b - i) - 1, M))
        if C[i] == -1:
            C[i] = []
        C[i].append(Temp)
    return C


def compensate(Compen, FFs):
    """Compensate the feedback functions in descending order
    Parameters
    ----------
    Compen : list
        compensation list
    FFs : list
        feedback functions of the Galois NLFSR
    Return
    ----------
    FFAll : list
        feedback functions after compensation
    """
    FFAll = copy.deepcopy(FFs)
    CompenStart = 0
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
            FFAll = sort_function(FFAll)

    return FFAll


def compensation_FF(FFGal, n):
    """Construct the feedback functions for the Fibonacci NLFSR
    Parameters
    ----------
    FFGal : list
        feedback functions of the Galois NLFSR
    n : int
        size of the NLFSR
    Return
    ----------
    FFAll : list
        feedback functions after compensation
    C_All : list
        compensation list
    """
    FFTemp = copy.deepcopy(FFGal)
    B = [] #The starting position of the shifted monomials
    MGal_All = [] #The monomials to be shifted
    MFib_All = [] #The monomials after shifting
    C_All = [-1] * n
    for i in range(n - 2, -1, -1): #Consider all the feedback functions except the last bit
        if len(FFTemp[i]) > 1:
            B.append(i)
            MGal = copy.deepcopy(FFTemp[i][1:])
            MFib = []
            Ci = [-1] * n
            for j in range(len(MGal)):
                Temp1 = list(map(lambda x: x - i - 1, MGal[j]))
                MFib.append(Temp1)
                Ci = compensation_list(i, MGal[j], Ci)
            MGal_All.append(MGal)
            MFib_All.append(MFib)
            for j in range(n): # Combine all the compensation sets
                if C_All[j] != -1 and Ci[j] != -1:
                    C_All[j] = C_All[j] + Ci[j]
                elif C_All[j] == -1:
                    C_All[j] = copy.deepcopy(Ci[j])
            FFTemp[i] = [[(i + 1) % n]] # Remove the shifted monomials gi in fi
            for j in range(n):
                if Ci[j] != -1:
                    Ci[j].insert(0, [j])
            FFTemp = compensate(Ci, FFTemp)
    for i in range(n):
        if C_All[i] != -1:
            C_All[i].insert(0, [i])
    FFAll = copy.deepcopy(FFTemp)
    for i in range(len(MFib_All)): #Add the shifted monomials to f_n-1
        FFAll[-1] = FFAll[-1] + MFib_All[i]
    return FFAll, C_All


def compensation_Fz(Z, CBackw):
    """Compensate the output function
    Parameters
    ----------
    Z : list
        output function of the Galois NLFSR
    CBackw : list
        compensation list
    Return
    ----------
    ZCompenBackw : list
        output function of the Fibonacci NLFSR
    """

    ZCompenBackw = []
    for i in range(len(Z)):
        Zi = copy.deepcopy(Z[i])  # Zi = [6, 4]
        Zi_temp = []
        for j in range(len(Zi)):
            if CBackw[Zi[j]] != -1:
                Zi_temp.append(CBackw[Zi[j]])
            else:
                Zi_temp.append([[Zi[j]]])
        Zi_temp = list(itertools.product(*Zi_temp))
        Zi_temp = [list(e) for e in Zi_temp]
        for j in range(len(Zi_temp)):
            Zi_temp[j] = [item for sublist in Zi_temp[j] for item in sublist]
            Zi_temp[j] = list(set(Zi_temp[j]))
        ZCompenBackw = ZCompenBackw + Zi_temp
    ZCompenBackw = sort_function([ZCompenBackw])
    ZCompenBackw = ZCompenBackw[0]

    return ZCompenBackw


def compute_IV(N0Origi, Compen, n):
    """Compute the initial state for the Galois NLFSR
    Parameters
    ----------
    N0Origi : list
        initial state of the Galois NLFSR, for example, [0, 1, 0, 1, 0, 0, 0]
    Compen : list
        compensation list, for example, [[[0], [1], [1], [2]], [[1], [2], [2], [3]], [[2], [3]], [[3], [4]], -1, -1, -1]
    n : int
        size of the NLFSR
    Return
    ----------
    N0Trans : list
        initial state of the Galois NLFSR, for example, [1, 0, 1, 1, 0, 0, 0]
    """

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



# ============================== An example below ==============================
# The size of the NLFSR
n = 7
# The number of rounds you want to run the NLFSR (The length of output sequence)
R = 1000

# The feedback functions of the Galois NLFSR
FFGal = []
for i in range(n):
    f = []
    f.insert(0, [(i + 1) % n])
    FFGal.append(f)
FFGal[1] = [[2], [3]]
FFGal[3] = [[4], [5]]
FFGal[6] = [[0], [4, 5]]
FFGal_origin = copy.deepcopy(FFGal)

# The output function of the Galois NLFSR
ZGal = [[2]]

# Generate a random initial state for the Galois NLFSR
N0Gal = []
for i in range(n):
    N0Gal.append(random.randint(0, 1))
# or set an initial state
N0Gal = [0, 1, 0, 1, 0, 0, 0]

OutGal, StatesGal = Gal_NLFSR(R, FFGal, ZGal, N0Gal)

# Compensate the feedback functions
FFFib, CompenAll = compensation_FF(FFGal, n)

# The final compensation list
CFinal = copy.deepcopy(CompenAll)
for i in range(n):
    if CompenAll[i] != -1:
        CFinal[i].append([i])
    else:
        CFinal[i] = 0
    CFinal = sort_function(CFinal)

# Construct output function and run the Fibonacci NLFSR
ZFib = compensation_Fz(ZGal, CompenAll)

# Compute the initial state for the Fibonacci NLFSR
N0Fib = compute_IV(N0Gal, CompenAll, n)

# Run the Fibonacci NLFSR and compare the output sequences
OutFib, StateFib = Fib_NLFSR(R, FFFib, ZFib, N0Fib)

if OutGal == OutFib:
    print('Transformation from Type-IV Galois NLFSR to Fibonacci NLFSR succeeds!')
    print('Feedback functions of the Galois NLFSR are: ', FFGal_origin)
    print('Output function of the Galois NLFSR is: ', ZGal)
    print('Initial state of the Galois NLFSR is: ', N0Gal)
    print('The combined compensation list is: ', CFinal)
    print('Feedback functions of the Fibonacci NLFSR are: ', FFFib)
    print('Output function of the Fibonacci NLFSR is: ', ZFib)
    print('Initial state of the Fibonacci NLFSR is: ', N0Fib)
    print('Output sequence is: ', OutFib)
else:
    print('Fail!')
