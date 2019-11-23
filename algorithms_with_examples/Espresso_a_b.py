"""The implementation of comparing Espresso, Espresso-a and Espresso-b
"""

import ast
import copy
from functools import reduce


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
    output = []
    for i in range(R):
        # First, compute the outputs
        Z = copy.deepcopy(ZF)
        for j in range(len(Z)):
            Z[j] = list(map(lambda x: N[x], Z[j]))

        for j in range(len(Z)):
            Z[j] = reduce(lambda p, q: p & q, Z[j])
        Z = reduce(lambda p, q: p ^ q, Z)
        output.append(Z)

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

    return output


# Then size of the Galois NLFSR
n = 256
# The number of output bits
R = 1000
# The initial state of the Espresso cipher
N0Gal= [0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1,
        0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1,
        1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1,
        1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0,
        1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1,
        0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
        1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
        0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1,
        1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0,
        0, 0, 0, 0, 1, 1]

# Generate output sequence of Espresso Cipher
FFGal = []
for i in range(n):
    f = []
    f.insert(0, [(i + 1) % n])
    FFGal.append(f)
FFGal[193] = [[194], [12, 121]]
FFGal[197] = [[198], [29, 52, 72, 99]]
FFGal[201] = [[202], [8, 103]]
FFGal[205] = [[206], [5, 80]]
FFGal[209] = [[210], [6, 64]]
FFGal[213] = [[214], [4, 45]]
FFGal[217] = [[218], [3, 32]]
FFGal[231] = [[232], [50, 159], [189]]
FFGal[235] = [[236], [67, 90, 110, 137]]
FFGal[239] = [[240], [46, 141], [117]]
FFGal[243] = [[244], [43, 118], [103]]
FFGal[247] = [[248], [44, 102], [40]]
FFGal[251] = [[252], [42, 83], [8]]
FFGal[255] = [[0], [41, 70]]
ZGal = [[80], [99], [137], [227], [222], [187], [243, 217], [247, 231], [213, 235], [255, 251],
        [181, 239], [174, 44], [164, 29], [255, 247, 243, 213, 181, 174]]
OutGal= Gal_NLFSR(R, FFGal, ZGal, N0Gal)

# Generate output sequence of the Espresso-a Cipher in reference [WL17]
N0Gal_a = [0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1,
           0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1,
           1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1,
           1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0,
           1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1,
           0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1,
           0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
           1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1,
           0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0,
           1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
           1, 1, 1, 0, 1, 0]
FFGal_a = []
for i in range(n):
    f = []
    f.insert(0, [(i + 1) % n])
    FFGal_a.append(f)
FFGal_a[254] = [[255], [40, 69]]
FFGal_a[255] = [[0], [12], [48], [115], [133], [213], [41, 70]]
ZGal_a = [[80], [99], [137], [227], [222], [187], [243, 217], [247, 231], [213, 235], [255, 251],
          [181, 239], [174, 44], [164, 29], [255, 247, 243, 213, 181, 174]]
OutGal_a = Gal_NLFSR(R, FFGal_a, ZGal_a, N0Gal_a)

# Generate output sequence of the Espresso-b Cipher transformed by using our proposed algorithm
# in Theorem 4 and Theorem 2
N0Gal_b = [0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1,
           0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1,
           1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1,
           1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0,
           1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1,
           0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1,
           0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
           1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1,
           0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0,
           1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
           1, 1, 1, 0, 1, 0]
FFGal_b = []
for i in range(n):
    f = []
    f.insert(0, [(i + 1) % n])
    FFGal_b.append(f)
FFGal_b[254] = [[255], [40, 69]]
FFGal_b[255] = [[0], [12], [48], [115], [133], [213], [41, 70]]

# read ZGal_b from file
with open('ZGal_b.txt') as f:
    ZGal_b = ast.literal_eval(f.read())

OutGal_b = Gal_NLFSR(R, FFGal_b, ZGal_b, N0Gal_b)

if OutGal == OutGal_a:
    print('Espresso is equivalent to Espresso-a!')
elif OutGal == OutGal_b:
    print('Espresso is equivalent to Espresso-b!')
