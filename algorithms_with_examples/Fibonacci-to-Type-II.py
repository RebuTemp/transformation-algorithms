"""
The implementation of the Fibonacci-to-Type-II transformation algorithm.
"""

import ast
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
    for i in range(len(FF)):
        tempx = FF[i]
        tempy = []
        for sublist in tempx:
            if sublist not in tempy:
                tempy.append(sublist)
            else:
                tempy.remove(sublist)
        FF_Sort.append(tempy)

    return FF_Sort


def compute_B(MFTG, n):
    """Compute possible positions to shift monomials
    Parameters
    ----------
    MFTG : list
        monomials to be shifted, for example, [[1], [2], [1, 3]]
    n : int
        size of the NLFSR
    Return
    ----------
    BFTG : list
        positions to shift monomials, for example, [3, 2, 3]
    """

    BFTG = []
    for i in range(len(MFTG)):
        m = copy.deepcopy(MFTG[i])
        if min(m) <= n - 1 - t:
            b = n - 1 - min(m)
        else:
            b = t
        BFTG.append(b)
    return BFTG


def compensation_list_FTG(MFTG, BFTG, n):
    """Construct compensation list
    Parameters
    ----------
    MFTG : list
        monomials to be shifted, for example, [[1], [2], [1, 3]]
    BFTG : list
        positions to shift monomials, for example, [3, 2, 3]
    n : int
        size of the NLFSR
    Return
    ----------
    CFTG : list
        compensation list, for example, [-1, -1, -1, [[1]]]
    """

    CFTG = [-1] * n
    for i in range(len(MFTG)):
        C = [-1] * n
        m = copy.deepcopy(MFTG[i])
        b = copy.deepcopy(BFTG[i])
        for j in range(b+1, n):
            CTemp = list(map(lambda x: x - n + j, m))
            C[j] = copy.deepcopy(CTemp)

        for j in range(n):
            if CFTG[j] == -1 and C[j] != -1:
                CFTG[j] = copy.deepcopy([C[j]])
            elif CFTG[j] != -1 and C[j] != -1:
                CFTG[j].append(C[j])

    return CFTG


def shift_FTG(FF, MFTG, BFTG, n):
    """Shift monomials
    Parameters
    ----------
    FF : list
        feedback functions of the Fibonacci NLFSR, for example, [[[1]], [[2]], [[3]], [[0], [1], [2], [1, 3]]]
    MFTG : list
        monomials to be shifted, for example, [[1], [2], [1, 3]]
    BFTG : list
        positions to shift monomials, for example, [3, 2, 3]
    n : int
        size of the NLFSR
    Return
    ----------
    FFShiftedFTG : list
        feedback functions after shifting monomials, for example, [[[1]], [[2]], [[3], [1]], [[0], [1], [1, 3]]]
    """

    FFShiftedFTG = copy.deepcopy(FF)
    for i in range(len(MFTG)):
        MTemp = list(map(lambda x: x - (n-1-BFTG[i]), MFTG[i]))
        FFShiftedFTG[n-1].append(MFTG[i])
        FFShiftedFTG[BFTG[i]].append(MTemp)
    FFShiftedFTG = sort_function(FFShiftedFTG)

    return FFShiftedFTG


def compensation_FTG(Z, CFTG, n):
    """Compensate the feedback functions and output function
    Parameters
    ----------
    FF : list
        feedback functions of the Fibonacci NLFSR, for example, [[[1]], [[2]], [[3]], [[0], [1], [2], [1, 3]]]
    Z : list
        output function of the Fibonacci NLFSR, for example, [[3]]
    CFTG : list
        compensation list, for example, [-1, -1, -1, [[1]]]
    n : int
        size of the NLFSR
    Return
    ----------
    FFCompenFTG : list
        feedback functions after compensation, for example, [[[1]], [[2]], [[3], [1]], [[0], [1, 3]]]
    ZCompenFTG : list
        output function after compensation, for example, [[1], [3]]
    """

    #FFCompenFTG = copy.deepcopy(FF)
    ZCompenFTG = copy.deepcopy(Z)
    #FFCompenFTG.append(Z)
    FFCompenFTG = copy.deepcopy([Z])
    for i in range(n-1, -1, -1):
        if CFTG[i] != -1:
            for j in range(1):
                FFj = copy.deepcopy(FFCompenFTG[j])
                for k in range(len(FFCompenFTG[j])):
                    if j < n and FFCompenFTG[j][k] == [(j + 1) % n]:
                        continue
                    FFjk = copy.deepcopy(FFCompenFTG[j][k])
                    if i in FFjk:
                        if FFjk in FFj:
                            FFj.remove(FFjk)
                        if len(FFCompenFTG[j][k]) == 1:
                            temp = copy.deepcopy(CFTG[i]+[[i]])
                        else:
                            FFjk = copy.deepcopy(FFCompenFTG[j][k])
                            FFjk.remove(i)
                            temp = list(map(lambda x: list(set(x+FFjk)), CFTG[i]+[[i]]))
                        FFj = FFj + temp
                FFCompenFTG[j] = copy.deepcopy(FFj)
            FFCompenFTG = sort_function(FFCompenFTG)
            ZCompenFTG = copy.deepcopy(FFCompenFTG[-1])

    return ZCompenFTG


def compute_IV_FTG(N0, CFTG, n):
    """Compute the initial state for the Galois NLFSR
    Parameters
    ----------
    N0 : list
        initial state of the Fibonacci NLFSR, for example, [1, 0, 1, 1]
    CFTG : list
        compensation list, for example, [-1, -1, -1, [[1]]]
    n : int
        size of the NLFSR
    Return
    ----------
    N0Gal : list
        initial state of the Galois NLFSR, for example, [1, 0, 1, 1]
    """
    N0Gal = copy.deepcopy(N0)
    CTemp = copy.deepcopy(CFTG)
    for i in range(n):
        if CTemp[i] != -1:
            gi = CTemp[i] + [[i]]
            for j in range(len(gi)):
                gi[j] = list(map(lambda x: N0[x], gi[j]))
            for j in range(len(gi)):
                gi[j] = reduce(lambda p, q: p & q, gi[j])
            gi = reduce(lambda p, q: p ^ q, gi)
            N0Gal[i] = gi

    return N0Gal


# ============================== An example below ==============================
# The size of the NLFSR
n = 4
# The number of rounds you want to run the NLFSR (The length of output sequence)
R = 1000

# The feedback function f_3 of the Fibonacci NLFSR
F = [[0], [1], [3], [1, 2]]
# All the feedback functions of the Fibonacci NLFSR
FF = []
for j in range(n - 1):
    f = [[j + 1]]
    FF.append(f)
FF.append(F)

# Output function of the Fibonacci NLFSR
Z = [[2, 3]]

# Generate a random initial state for the Fibonacci NLFSR
N0 = []
for i in range(n):
    N0.append(random.randint(0, 1))
N0 = [0, 0, 0, 1]

print('FFFib = ', FF)
print('ZFib = ', Z)
print('N0Fib = ', N0)
print('\n')

# Monomials to be shifted (except the monomial x_0)
M = copy.deepcopy(F[1:])
# Compute possible positions to shift monomials
T = []
for i in range(len(F)):
    if len(F[i]) > 0:
        T.append(max(F[i]) - min(F[i]))
t = max(T)
B = compute_B(M, n)

# All the combination of possible positions
BAll = []
for i in range(len(B)):
    BTemp = list(range(n-1, B[i]-1, -1))  # n-1
    BAll.append(BTemp)
BAll = list(itertools.product(*BAll))
BAll = [list(e) for e in BAll]

# For each combination of possible positions, we transform the Fibonacci NLFSR to a Galois NLFSR
for i in range(len(BAll)):
    FFFib = copy.deepcopy(FF)
    ZFib = copy.deepcopy(Z)
    N0Fib = copy.deepcopy(N0)
    MFTG = copy.deepcopy(M)
    BFTG = copy.deepcopy(BAll[i])

    # First, we shift monomials
    FFGal = shift_FTG(FFFib, MFTG, BFTG, n)

    # Second, we construct the compensation list
    CFTG = compensation_list_FTG(MFTG, BFTG, n)

    # Third, we compensate the feedback functions and output function by using the compensation list
    ZGal = compensation_FTG(ZFib, CFTG, n)

    # Last, we compute the initial state for the transformed Galois NLFSR
    N0Gal = compute_IV_FTG(N0, CFTG, n)

    OutFib, StateFib = Fib_NLFSR(R, FFFib, ZFib, N0Fib)
    OutGal, StateGal = Gal_NLFSR(R, FFGal, ZGal, N0Gal)
    if OutFib == OutGal:
        print('===== Possible Transformation: %s' % (i + 1))
        print('Transformation from Fibonacci NLFSR to Galois NLFSR succeeds!')
        print('We shift monomials to feedback functions with indexes: ', BFTG)
        print('The compensation list is: ', CFTG)
        print('Feedback functions of the Galois NLFSR are: ', FFGal)
        print('Output function of the Galois NLFSR is: ', ZGal)
        print('Initial state of the Galois NLFSR is: ', N0Gal)
        print('------------------------------------------------------------\n')
    else:
        print('Transformation fails because the Galois NLFSR is not uniform!')
        print('Feedback functions of the Galois NLFSR are: ', FFGal)
        print('------------------------------------------------------------\n')

