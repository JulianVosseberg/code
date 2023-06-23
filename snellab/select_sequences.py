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
if len(queries) == 0 or queries == ['']:
    sys.exit("Error: no queries found.")

found = []

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
                        found.append(query)
                        if args.t:
                            queries.remove(query)
                        break
            else:
                if to_include:
                    print(line)
    else:
        sequences_dict = {}
        for line in fasta_file:
            if line[0] == '>':
                seqid = line[1:].rstrip()
                if seqid in sequences_dict:
                    sys.exit(f'Error: multiple headers with the same name ({seqid}). Make all sequence IDs unique.')
                sequences_dict[seqid] = ''
            else:
                sequences_dict[seqid] += line
        for seqid, sequence in sequences_dict.items():
            to_include = True
            if args.f:
                if seqid in queries:
                    to_include = False
                    if args.t:
                        query = seqid
            else:
                for query in queries:
                    if seqid.startswith(query):
                        to_include = False
                        break
            if to_include:
                print(f'>{seqid}\n{sequence}', end = '')
            else:
                found.append(query)
                if args.t:
                    queries.remove(query)
if args.t and len(queries) > 0:
    sys.stderr.write(f'Warning: the following queries were not found: {",".join(queries)}\n')
elif not args.t:
    not_found = set(queries).difference(set(found))
    if len(not_found) > 0:
        sys.stderr.write(f'Warning: the following queries were not found: {",".join(not_found)}\n')
