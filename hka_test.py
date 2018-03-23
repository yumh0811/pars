#/usr/bin/env python3

# python hka_test.py 9 8 210 18 414 79 4052 324 --name Adh --size 81

# Mathematica
# Solve[{0 == i * (e * y + f * z) - (a + b), 0 == (g * y + h * z ) * (x + 1) - (c + d), 0 == g * y * (x + 1) + i * e * y - (c + a)}, {x, y, z}]

import argparse
from sympy import *
from sympy.core.cache import *
import itertools

parser = argparse.ArgumentParser(description='Original HKA test')
parser.add_argument(
    'diff',
    type=int,
    nargs=4,
    help='four integers of differences: S1A S2A D1 D2')
parser.add_argument(
    'comp',
    type=int,
    nargs=4,
    help='four integers of comparable: CS1A CS2A CD1 CD2')
parser.add_argument('--name', default='locus', help='the name of test locus')
parser.add_argument('--size', type=int, default=81, help='population size')

args = parser.parse_args()

# clear_cache()

# constant

# var('S1A S2A D1 D2')
S1A = args.diff[0]
S2A = args.diff[1]
D1 = args.diff[2]
D2 = args.diff[3]

# var('CS1A CS2A CD1 CD2')
CS1A = args.comp[0]
CS2A = args.comp[1]
CD1 = args.comp[2]
CD2 = args.comp[3]

nA = args.size
f = 1

# var('CnA')
var('j')
CnA = Sum(1 / j, (j, 1, nA - 1)).doit()
CnA2 = Sum(1 / j**2, (j, 1, nA - 1)).doit()

# var
That, theta1, theta2, X2 = var('That theta1 theta2 X2')

ES1A = CS1A * theta1 * CnA
ES2A = CS2A * theta2 * CnA
ED1 = CD1 * theta1 * (That + (1 + f) / 2)
ED2 = CD2 * theta2 * (That + (1 + f) / 2)

# equation system
E1 = CnA * (CS1A * theta1 + CS2A * theta2) - (S1A + S2A)
E2 = (CD1 * theta1 + CD2 * theta2) * (That + 1) - (D1 + D2)
E3 = CD1 * theta1 * (That + 1) + CnA * CS1A * theta1 - (D1 + S1A)
E4 = (S1A - ES1A)**2 / (ES1A + (CS1A * theta1)**2 * CnA2) + \
     (S2A - ES2A)**2 / (ES2A + (CS2A * theta2)**2 * CnA2) + \
     (D1  - ED1)**2  / (ED1  + (CD1 * theta1*(1+f)/2)**2) + \
     (D2  - ED2)**2  / (ED2  + (CD2 * theta2*(1+f)/2)**2) - \
     X2

sols = solve([E1, E2, E3, E4], [That, theta1, theta2, X2])

flatten = list(itertools.chain.from_iterable(sols))
out = "\t".join(map(str, flatten))
print("%s\t%s\t%s\n" % (args.name, args.size, out) )
