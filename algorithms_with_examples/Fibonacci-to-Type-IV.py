"""
The implementation of the Fibonacci-to-Type-IV transformation algorithm.
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

def Compu_B_Backw(MBackw, n):
    """Compute possible positions to shift monomials
    Parameters
    ----------
    MBackw : list
        monomials to be shifted, for example, [[1], [2], [1]]
    n : int
        size of the NLFSR
    Return
    ----------
    BBackw : list
        positions to shift monomials, for example, [1, 1, 3]
    """

    BBackw = []
    for i in range(len(MBackw)):
        m = copy.deepcopy(MBackw[i])
        b = n - 1 - max(m) - 1
        BBackw.append(b)
    return BBackw

def Compen_List_Backw(MBackw, BBackw, n):
    """Construct compensation list
    Parameters
    ----------
    MBackw : list
        monomials to be shifted, for example, [[1], [2], [1]]
    BBackw : list
        positions to shift monomials, for example, [1, 1, 3]
    n : int
        size of the NLFSR
    Return
    ----------
    CBackw : list
        compensation list, for example, [[[1], [2], [1]], [[2], [3], [2]], [[3]], [[4]], -1, -1, -1]
    """

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
    """Shift monomials
    Parameters
    ----------
    FF : list
        feedback functions of the Fibonacci NLFSR, for example, [[[1]], [[2]], [[3]], [[4]], [[5]], [[6]], [[0], [4, 5], [1], [2], [1]]]
    MBackw : list
        monomials to be shifted, for example, [[1], [2], [1]]
    BBackw : list
        positions to shift monomials, for example, [1, 1, 3]
    n : int
        size of the NLFSR
    Return
    ----------
    FFShiftedBackw : list
        feedback functions after shifting monomials, for example, [[[1]], [[2], [3], [4]], [[3]], [[4], [5]], [[5]], [[6]], [[0], [4, 5]]]
    """

    FFShitedBackw = copy.deepcopy(FF)
    for i in range(len(MBackw)):
        MTemp = list(map(lambda x: x + BBackw[i] + 1, MBackw[i]))
        FFShitedBackw[n-1].append(MBackw[i])
        FFShitedBackw[BBackw[i]].append(MTemp)
    FFShitedBackw = sort_function(FFShitedBackw)
    return FFShitedBackw

def Compensation_Backw(FF, Z, CBackw, n):
    """Compensate the feedback functions and output function
    Parameters
    ----------
    FF : list
        feedback functions of the Fibonacci NLFSR, for example, [[[1]], [[2], [3], [4]], [[3]], [[4], [5]], [[5]], [[6]], [[0], [4, 5]]]
    Z : list
        output function of the Fibonacci NLFSR, for example, [[2], [3]]
    CBackw : list
        compensation list, for example, [[[1], [2], [1]], [[2], [3], [2]], [[3]], [[4]], -1, -1, -1]
    n : int
        size of the NLFSR
    Return
    ----------
    FFCompenBackw[:-1] : list
        feedback functions after compensation, for example, [[[1]], [[2], [3]], [[3]], [[4], [5]], [[5]], [[6]], [[0], [4, 5]]]
    ZCompenBackw : list
        output function after compensation, for example, [[2]]
    """

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
            FFCompenBackw = sort_function(FFCompenBackw)
            ZCompenBackw = copy.deepcopy(FFCompenBackw[-1])
    return FFCompenBackw[:-1], ZCompenBackw

def Compute_IV(N0, C, n):
    """Compute the initial state for the Galois NLFSR
    Parameters
    ----------
    N0 : list
        initial state of the Fibonacci NLFSR, for example, [1, 0, 1, 1, 0, 0, 0]
    C : list
        compensation list, for example, [[[1], [2], [1]], [[2], [3], [2]], [[3]], [[4]], -1, -1, -1]
    n : int
        size of the NLFSR
    Return
    ----------
    N0Gal : list
        initial state of the Galois NLFSR, for example, [0, 1, 0, 1, 0, 0, 0]
    """

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


# ============================== An example below ==============================
# The size of the NLFSR
n = 7
# The number of rounds you want to run the NLFSR (The length of output sequence)
R = 1000

# The feedback function f_{n-1} of the Fibonacci NLFSR
F = [[0], [4, 5], [1], [2], [1]]
# All the feedback functions of the Fibonacci NLFSR
FF = []
for j in range(n - 1):
    f = [[j + 1]]
    FF.append(f)
FF.append(F)

# Output function of the Fibonacci NLFSR
Z = [[2], [3]]

# Generate a random initial state for the Fibonacci NLFSR
N0 = []
for i in range(n):
    N0.append(random.randint(0, 1))
# or set an initial state
N0 = [1, 0, 1, 1, 0, 0, 0]

print('FFFib = ', FF)
print('ZFib = ', Z)
print('N0Fib = ', N0)
print('\n')

# Choose the monomials to be shifted
M = copy.deepcopy(F[2:])

# Calculate all the combinations of possible positions
B = Compu_B_Backw(M, n)
BAll = []
for i in range(len(B)):
    BTemp = list(range(0, B[i]+1, 1)) #n-1
    BAll.append(BTemp)
BAll = list(itertools.product(*BAll))
BAll = [list(e) for e in BAll]

# For each combination of possible positions, we transform the Fibonacci NLFSR to a Type-IV Galois NLFSR
for i in range(len(BAll)):
    FFFib = copy.deepcopy(FF)
    ZFib = copy.deepcopy(Z)
    N0Fib = copy.deepcopy(N0)
    MBackw = copy.deepcopy(M)
    BBackw = copy.deepcopy(BAll[i])

    # Construct the compensation list
    CBackw = Compen_List_Backw(MBackw, BBackw, n)

    # Shift the monomials
    FFShiftedBackw = Shifting_Backw(FFFib, MBackw, BBackw, n)

    # Compensate the feedback functions and the output function
    FFGal, ZGal = Compensation_Backw(FFShiftedBackw, ZFib, CBackw, n)

    # Compute the initial state
    N0Gal = Compute_IV(N0, CBackw, n)

    # Compare the output sequences
    OutFib, StateFib = Fib_NLFSR(R, FFFib, ZFib, N0Fib)
    OutGal, StateGal = Gal_NLFSR(R, FFGal, ZGal, N0Gal)
    if OutFib == OutGal:
        print('===== Possible Transformation: %s' % (i + 1))
        print('Transformation from Fibonacci NLFSR to Galois NLFSR succeeds!')
        print('We shift monomials to feedback functions with indexes: ', BBackw)
        print('The compensation list is: ', CBackw)
        print('Feedback functions of the Galois NLFSR are: ', FFGal)
        print('Output function of the Galois NLFSR is: ', ZGal)
        print('Initial state of the Galois NLFSR is: ', N0Gal)
        print('Output sequence is: ', OutGal)
        print('------------------------------------------------------------\n')
    else:
        print('Transformation fails because the Galois NLFSR is not uniform!')
        print('Feedback functions of the Galois NLFSR are: ', FFGal)
        print('------------------------------------------------------------\n')








