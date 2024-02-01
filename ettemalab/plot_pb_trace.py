#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 15:49:57 2024

@author: julian
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style = 'whitegrid')

parser = argparse.ArgumentParser(description = "This script reads the trace files of the PhyloBayes runs and dusplays the sampled parameter values.")
parser.add_argument("prefix", help = "prefix (may include path) of the chains")
parser.add_argument("-c", metavar = "chains", help = "number of chains (default = 2)", type = int, default = 2)
parser.add_argument("-x", metavar = "xmax", help = "maximum value for x axis (number of cycles) to display", type = int)
parser.add_argument("-b", metavar = "burnin", help = "minimum value for x axis (number of cycles) to display", type = int, default = 0)
parser.add_argument("-o", metavar = "outdir", help = "output directory for plot (current directory by default)", default = ".")
args = parser.parse_args()

for i in range(1, args.c + 1):
    if args.b:
        df = pd.read_table(f'{args.prefix}{i}.trace', skiprows = range(1, args.b))
    else:
        df = pd.read_table(f'{args.prefix}{i}.trace')
    df['chain'] = f'chain{i}'
    if i == 1:
        tracedf = df
    else:
        tracedf = pd.concat([tracedf, df], ignore_index = True)

ncol_to_plot = tracedf.shape[1] - 4
cols_to_plot = list(tracedf.columns[3:-1])
if ncol_to_plot == 8:
    nrows = 4
    ncols = 2
    figsize = (12, 8)
elif ncol_to_plot == 10:
    nrows = 5
    ncols = 2
    figsize = (12, 10)
else:
    sys.exit("Error: add a new way to plot this number of parameters!")
fig, axs = plt.subplots(nrows, ncols, figsize = figsize)
k = 0
for i in range(nrows):
    for j in range(ncols):
        column = tracedf.columns[k + 3]
        legend = False
        if k == 0:
            legend = 'full'
        sns.lineplot(data = tracedf,
                     x = '#cycle',
                     y = column,
                     hue = 'chain',
                     legend = legend,
                     ax = axs[i,j])
        if args.x:
            axs[i,j].set_xlim(right = args.x)
#        else:
#            axs[i,j].set_xlim(args.b)
        k += 1
sns.despine(bottom = True, left = True)
if args.x and args.b:
    outname = f'traceplots_xmin_{args.b}_xmax_{args.x}.pdf'
elif args.b:
    outname = f'traceplots_xmin_{args.b}.pdf'
elif args.x:
    outname = f'traceplots_xmax_{args.x}.pdf'
else:
    outname = 'traceplots.pdf'
fig.savefig(f'{args.o}/{outname}', bbox_inches = 'tight')