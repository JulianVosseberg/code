#! /usr/bin/env python3

import sys
from random import sample

def usage():
    sys.exit("\tUsage: get_random_selection.py <#> <item1> <item2> ...\n\
    \tThis script selects a specified number of items from the given arguments.")

if len(sys.argv) < 4:
    usage()

number = int(sys.argv[1])
population = sys.argv[2:]

if len(population) <= number:
    print('Population size is not larger than the specified number!\n', file = sys.stderr)
    usage()

selection = sample(population, number)
print(' '.join(selection))
