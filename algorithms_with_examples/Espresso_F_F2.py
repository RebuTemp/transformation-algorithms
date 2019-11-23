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
N0Gal_G= [1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1]

# Generate output sequence of Espresso Cipher
FFGal_G = []
for i in range(n):
    f = []
    f.insert(0, [(i + 1) % n])
    FFGal_G.append(f)
FFGal_G[193] = [[194], [12, 121]]
FFGal_G[197] = [[198], [29, 52, 72, 99]]
FFGal_G[201] = [[202], [8, 103]]
FFGal_G[205] = [[206], [5, 80]]
FFGal_G[209] = [[210], [6, 64]]
FFGal_G[213] = [[214], [4, 45]]
FFGal_G[217] = [[218], [3, 32]]
FFGal_G[231] = [[232], [50, 159], [189]]
FFGal_G[235] = [[236], [67, 90, 110, 137]]
FFGal_G[239] = [[240], [46, 141], [117]]
FFGal_G[243] = [[244], [43, 118], [103]]
FFGal_G[247] = [[248], [44, 102], [40]]
FFGal_G[251] = [[252], [42, 83], [8]]
FFGal_G[255] = [[0], [41, 70]]
ZGal_G = [[80], [99], [137], [227], [222], [187], [243, 217], [247, 231], [213, 235], [255, 251],
        [181, 239], [174, 44], [164, 29], [255, 247, 243, 213, 181, 174]]
OutGal_G = Gal_NLFSR(R, FFGal_G, ZGal_G, N0Gal_G)


# Generate output sequence of the Galois NLFSR F in reference [DH17]
N0Gal_F = [1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0]

FFGal_F = []
for i in range(n):
    f = []
    f.insert(0, [(i + 1) % n])
    FFGal_F.append(f)
FFGal_F[217] = [[218], [3, 32], [8, 49], [14, 72], [17, 92], [24, 119], [36, 145], [49, 72, 92, 119]]
FFGal_F[255] = [[0], [12], [48], [115], [133], [213], [41, 70], [46, 87], [52, 110], [55, 130], [62, 157], [74, 183], [87, 110, 130, 157]]
ZGal_F = [[80], [99], [137], [227], [222], [187], [243, 217], [247, 231], [213, 235], [255, 251],
          [181, 239], [174, 44], [164, 29], [255, 247, 243, 213, 181, 174]]
OutGal_F = Gal_NLFSR(R, FFGal_F, ZGal_F, N0Gal_F)

# Generate output sequence of the Galois NLFSR F2 transformed by using our proposed algorithm
N0Gal_F2 = [1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0]

FFGal_F2 = []
for i in range(n):
    f = []
    f.insert(0, [(i + 1) % n])
    FFGal_F2.append(f)
FFGal_F2[217] = [[218], [3, 32], [8, 49], [14, 72], [17, 92], [24, 119], [36, 145], [49, 72, 92, 119]]
FFGal_F2[255] = [[0], [12], [48], [115], [133], [213], [41, 70], [46, 87], [52, 110], [55, 130], [62, 157], [74, 183], [87, 110, 130, 157]]

# read ZGal_F2 from file
with open('ZGal_F2.txt') as f2:
    ZGal_F2 = ast.literal_eval(f2.read())

OutGal_F2 = Gal_NLFSR(R, FFGal_F2, ZGal_F2, N0Gal_F2)

if OutGal_G == OutGal_F:
    print('Espresso is equivalent to Espresso based on Galois NLFSR F!')
elif OutGal_G == OutGal_F2:
    print('Espresso is equivalent to Espresso based on Galois NLFSR F2!')