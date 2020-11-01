"""
The implementation of the Type-II-to-Fibonacci transformation algorithm.
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


def compensation_list_GTF(FF, n):
    """Construct compensation list
    Parameters
    ----------
    FF : list
        feedback functions of the Galois NLFSR
    n : int
        size of the NLFSR
    Return
    ----------
    Compen : list
        compensation list
    MGal : list
        Monomial after shifting
    B : list
        position where the monomial is shifted from
    """

    FFTemp = copy.deepcopy(FF)
    Compen = [-1] * n
    MGal = []
    B = []
    for i in range(n - 1):
        CompenTemp = [-1] * n
        if len(FFTemp[i]) > 1:
            for j in range(1, len(FFTemp[i])):
                MGal.append(FFTemp[i][j])
                B.append(i)
            CompenStart = i + 1
            FTemp = copy.deepcopy(FFTemp[i])
            CompenTemp[i + 1:] = list(map(lambda x: [], CompenTemp[i + 1:]))
            CompenTemp[i + 1] = copy.deepcopy(FTemp)
            for j in range(len(FTemp)):
                pTemp = copy.deepcopy(FTemp[j])
                for k in range(CompenStart + 1, n):
                    pTemp = list(map(lambda x: (x + 1) % n, pTemp))
                    CompenTemp[k].append(pTemp)
            break
    for j in range(len(CompenTemp)):
        if CompenTemp[j] != -1:
            if Compen[j] == -1:
                Compen[j] = copy.deepcopy(CompenTemp[j])
            else:
                Compen[j] = Compen[j] + CompenTemp[j][1:]

    return Compen, MGal, B


def compensation_GTF(FF, CompenAll, CompenPosition, n):
    """Construct the feedback functions for the Fibonacci NLFSR,
    only shifting process is needed so we comment out the compensation process
    Parameters
    ----------
    FF : list
        feedback functions of the Galois NLFSR
    CompenAll : list
        compensation list
    CompenPosition : list
        start position to compensate the feedback function
    n : int
        size of the NLFSR
    Return
    ----------
    FFAll : list
        feedback functions after compensation
    """

    FFAll = copy.deepcopy(FF)
    i = CompenPosition + 1
    for j in range(n):
        FFj = copy.deepcopy(FFAll[j])
        for k in range(len(FFAll[j])):
            if FFAll[j][k] == [(j + 1) % n]:
                continue
            FFjk = copy.deepcopy(FFAll[j][k])
            if i in FFjk:
                if FFjk in FFj:
                    FFj.remove(FFjk)
                # if len(FFAll[j][k]) == 1:
                    # temp = copy.deepcopy(CompenAll[i])
                # else:
                    # FFjk = copy.deepcopy(FFAll[j][k])
                    # FFjk.remove(i)
                    # temp = list(map(lambda x: list(set(x+FFjk)), CompenAll[i]))
                # FFj = FFj + temp
        FFAll[j] = copy.deepcopy(FFj)
    FFAll = sort_function(FFAll)

    return FFAll


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


def compute_IV_GTF(N0, Compen, n):
    """Compute the initial state for the Galois NLFSR
    Parameters
    ----------
    N0 : list
        initial state of the Galois NLFSR, for example, [0, 1, 0, 0]
    CFTG : list
        compensation list, for example, [-1, -1, [[2], [1]], [[3], [2], [0, 1]]]
    n : int
        size of the NLFSR
    Return
    ----------
    N0GalComp : list
        initial state of the Fibonacci NLFSR, for example, [0, 1, 1, 1]
    """

    N0GalComp = copy.deepcopy(N0)
    CompenTemp = copy.deepcopy(Compen)
    for i in range(n):
        if CompenTemp[i] != -1:
            gi = CompenTemp[i]
            for j in range(len(gi)):
                gi[j] = list(map(lambda x: N0GalComp[x], gi[j]))
            for j in range(len(gi)):
                gi[j] = reduce(lambda p, q: p & q, gi[j])
            gi = reduce(lambda p, q: p ^ q, gi)
            N0GalComp[i] = gi

    return N0GalComp


# ============================== An example below ==============================
# The size of the NLFSR
n = 4
# The number of rounds you want to run the NLFSR (The length of output sequence)
R = 1000

# The feedback functions of the Galois NLFSR
FFGal = []
for i in range(n):
    f = []
    f.insert(0, [(i + 1) % n])
    FFGal.append(f)
FFGal[1] = [[2], [1]]
FFGal[2] = [[3], [0, 1]]
FFGal[3] = [[0], [1]]
FFGal_origin = copy.deepcopy(FFGal)

# The output function of the Galois NLFSR
ZGal = [[1], [2], [0, 1], [0, 1, 2], [1, 3], [2, 3]]

# Generate a random initial state for the Galois NLFSR
N0Gal = []
for i in range(n):
    N0Gal.append(random.randint(0, 1))
# or set an initial state
N0Gal = [0, 1, 0, 0]

OutGal, StatesGal = Gal_NLFSR(R, FFGal, ZGal, N0Gal)

# Step 1: construct compensation list and feedback functions for the Fibonacci NLFSR
CompenAll = [-1] * n
MFib = []
MGalAll = []
BAll = []
for x in range(n - 1):
    Compen, MGal, B = compensation_list_GTF(FFGal, n)
    for i in range(len(B)):
        MGali = copy.deepcopy(MGal[i])
        MGali = list(map(lambda x: x + n - 1 - B[i], MGali))
        MFib.append(MGali)

    for i in range(n):
        if Compen[i] != -1:
            if CompenAll[i] == -1:
                CompenAll[i] = copy.deepcopy(Compen[i])
            else:
                CompenAll[i] = CompenAll[i] + Compen[i][1:]
    MGalAll = MGalAll + MGal
    BAll = BAll + B
    if B == []:
        break
    CompenPosition = B[0]
    FFGal[CompenPosition] = [x for x in FFGal[CompenPosition] if x not in MGal]
    FFGal = compensation_GTF(FFGal, CompenAll, CompenPosition, n)
    for i in range(B[0] + 1, n):
        if len(FFGal[i]) == 1:
            CompenPosition = i
            FFGal = compensation_GTF(FFGal, CompenAll, CompenPosition, n)
        else:
            break
FFGal[-1] = FFGal[-1] + MFib
FFGal = sort_function(FFGal)
FFFib = copy.deepcopy(FFGal)

# The final compensation list
CFinal = copy.deepcopy(CompenAll)
for i in range(n):
    if CompenAll[i] != -1:
        CFinal[i].append([i])
    else:
        CFinal[i] = 0
    CFinal = sort_function(CFinal)

# Step 2: Compensate the output function by using compensation list
ZFib = compensation_Fz(ZGal, CompenAll)

# Step 3: Compute the initial state for the Fibonacci NLFSR
N0Fib = compute_IV_GTF(N0Gal, CompenAll, n)

OutFib, StatesFib = Fib_NLFSR(R, FFFib, ZFib, N0Fib)

if OutGal == OutFib:
    print('Transformation from Type-II Galois NLFSR to Fibonacci NLFSR succeeds!')
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
