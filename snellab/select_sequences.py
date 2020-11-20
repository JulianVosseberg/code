#!/usr/bin/env python

# Load modules
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser(description = "This script extracts or removes a sequence or sequences of interest from a FASTA file.")
parser.add_argument("fasta", help = "FASTA file")
group = parser.add_mutually_exclusive_group()
group.add_argument('-q', metavar = "queries", help = "comma-separated queries (query1,query2,query3,...)")
group.add_argument('-l', metavar = "list", help = "file containing list of queries")
parser.add_argument('-m', help = "mode, default is to extract", choices = ('extract', 'remove'), default = 'extract')
parser.add_argument('-f', help = 'full match required', action = 'store_true')
parser.add_argument('-t', help = 'top hit (stop looking after first hit)', action = 'store_true')
args = parser.parse_args()

if args.q:
    queries = args.q.split(',')
elif args.l:
    with open(args.l, 'r') as query_list:
        queries = query_list.read().rstrip().split('\n')
if len(queries) == 0:
    sys.exit("Error: no queries found.")

if args.m == "remove":
    sys.exit("Sorry, remove mode not implemented yet...")

with open(args.fasta, 'r') as fasta_file:
    if args.m == "extract":
        to_include = False
        for line in fasta_file:
            line = line.rstrip()
            if line[0] == '>':
                if to_include:
                    to_include = False
                    if args.t and len(queries) == 0:
                        break
                for query in queries:
                    if args.f:
                        if line[1:] == query:
                            to_include = True
                    else:
                        if line[1:].startswith(query):
                            to_include = True
                    if to_include:
                        print(line)
                        if args.t:
                            queries.remove(query)
                        break
            else:
                if to_include:
                    print(line)
    else:
        to_include = True
        for line in fasta_file:
            line = line.rstrip()
            if line[0] == '>':
                for query in queries:
                    pass
            else:
                if to_include:
                    print(line)
