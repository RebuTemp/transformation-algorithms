"""
The implementation of the Type-I-to-Fibonacci transformation algorithm.
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


def compensation_list(FF, n):
    """Construct the feedback functions for the Fibonacci NLFSR
    Parameters
    ----------
    FF : list
        feedback functions of the Galois NLFSR
    n : int
        size of the NLFSR
    Return
    ----------
    FFAll : list
        feedback functions after compensation
    CompenAll : list
        combined compensation list
    """

    FFAll = copy.deepcopy(FF)
    CompenAll = [-1]*n
    for i in range(n-1):
        Compen = [-1]*n
        if len(FF[i]) > 1:
            FFAll[i] = [[i+1]]
            for j in range(i+1):
                Compen[j] = [n-1-i+j]
            for i in range(n):
                if Compen[i] != -1:
                    if CompenAll[i] == -1:
                        CompenAll[i] = [Compen[i]]
                    else:
                        CompenAll[i].append(Compen[i])
    CompenAll = sort_function(CompenAll)

    return FFAll, CompenAll


def compensation_GTF(F, CompenList):
    """Compensate the feedback function and the output function for the Fibonacci NLFSR
    Parameters
    ----------
    F : list
        a feedback function of the Galois NLFSR
    CompenAll : list
        compensation list
    n : int
        size of the NLFSR
    Return
    ----------
    FCompen : list
        the feedback function after compensation
    """
    C = copy.deepcopy(CompenList)
    FCompen = []
    for i in range(len(F)):
        Fi = copy.deepcopy(F[i])  # Fi = [0, 4]
        Fi_temp = []
        for j in range(len(Fi)):
            if C[Fi[j]] != -1:
                Fi_temp.append(C[Fi[j]])
                Fi_temp[j].append([Fi[j]])
            else:
                Fi_temp.append([[Fi[j]]])
        Fi_temp = list(itertools.product(*Fi_temp))
        Fi_temp = [list(e) for e in Fi_temp]
        for j in range(len(Fi_temp)):
            Fi_temp[j] = [item for sublist in Fi_temp[j] for item in sublist]
            Fi_temp[j] = list(set(Fi_temp[j]))
        FCompen = FCompen + Fi_temp
    FCompen = sort_function([FCompen])
    FCompen = FCompen[0]

    return FCompen


def compute_IV_GTF(N0, Compen, n):
    """Compute the initial state for the Galois NLFSR
    Parameters
    ----------
    N0 : list
        initial state of the Galois NLFSR, for example, [0, 1, 1, 0, 1]
    CFTG : list
        compensation list, for example, [[[4], [3], [1]], [[4], [2]], [[3]], [[4]], -1]
    n : int
        size of the NLFSR
    Return
    ----------
    N0Gal : list
        initial state of the Fibonacci NLFSR, for example, [0, 0, 0, 1, 1]
    """

    N0GalComp = copy.deepcopy(N0)
    CompenTemp = copy.deepcopy(Compen)
    for i in range(n-1,-1,-1):
        if CompenTemp[i] != -1:
            gi = CompenTemp[i]
            for j in range(len(gi)):
                gi[j] = list(map(lambda x: N0GalComp[x], gi[j]))
            for j in range(len(gi)):
                gi[j] = reduce(lambda p, q: p & q, gi[j])
            gi.append(N0[i])
            gi = reduce(lambda p, q: p ^ q, gi)
            N0GalComp[i] = gi

    return N0GalComp


# ============================== An example below ==============================
# The size of the NLFSR
n = 5
# The number of rounds you want to run the NLFSR (The length of output sequence)
R = 1000

# The feedback functions of the Galois NLFSR
FFGal = []
for i in range(n):
    f = []
    f.insert(0, [(i + 1) % n])
    FFGal.append(f)
FFGal[0] = [[1], [0], [1, 2]]
FFGal[1] = [[2], [0], [1, 2]]
FFGal[3] = [[4], [0], [1, 2]]
FFGal[4] = [[0], [1, 2]]
FFGal_origin = copy.deepcopy(FFGal)

# The output function of the Galois NLFSR
ZGal = [[2], [0, 4]]

# Generate a random initial state for the Galois NLFSR
N0Gal = []
for i in range(n):
    N0Gal.append(random.randint(0, 1))
# or set an initial state
N0Gal = [0, 1, 1, 0, 1]

OutGal, StatesGal = Gal_NLFSR(R, FFGal, ZGal, N0Gal)

# Step 1 and Step 2: Remove monomials and construct the combined compensation list
FFFib, CompenList = compensation_list(FFGal,n)

# Step 3 : Compensate the feedback function
FFFib[n-1] = compensation_GTF(FFFib[-1], CompenList)

# Step 4: Compensate the output function
ZFib = compensation_GTF(ZGal, CompenList)

# Step 5: Compute IV for the Fibonacci NLFSR
N0Fib = compute_IV_GTF(N0Gal,CompenList,n)

OutFib, StatesFib = Fib_NLFSR(R, FFFib, ZFib, N0Fib)

if OutGal == OutFib:
    print('Transformation from Type-I Galois NLFSR to Fibonacci NLFSR succeeds!')
    print('Feedback functions of the Galois NLFSR are: ', FFGal_origin)
    print('Output function of the Galois NLFSR is: ', ZGal)
    print('Initial state of the Galois NLFSR is: ', N0Gal)
    print('The combined compensation list is: ', CompenList)
    print('Feedback functions of the Fibonacci NLFSR are: ', FFFib)
    print('Output function of the Fibonacci NLFSR is: ', ZFib)
    print('Initial state of the Fibonacci NLFSR is: ', N0Fib)
    print('Output sequence is: ', OutFib)
else:
    print('Fail!')

