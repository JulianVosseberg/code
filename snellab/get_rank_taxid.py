#!/usr/bin/python

import sys; import getopt
from Bio import Entrez

def usage():
    print "\tUsage: get_rank_taxid.py [-h] -r <rank1,rank2,..> -l <taxids_list> -q <taxid1,taxid2,..>"
    print "\tThis script retrieves the desired rank for each NCBI taxonomy ID"
    print "\nUse -l or -q"

def get_rank(taxid, ranks):
    handle = Entrez.efetch(db = "Taxonomy", id = taxid, retmode = "xml")
    records = Entrez.read(handle)
    classification = []
    for rank in records[0]['LineageEx']:
        if rank["Rank"] in ranks:
            classification.append(rank["ScientificName"])
    return(classification)

optlist, args = getopt.getopt(sys.argv[1:], 'r:l:q:h')
opts = {}
ids_list = False; queries = False
for k,v in optlist:
    if k == "-h":
        usage(); sys.exit()
    else:
        opts[k] = v
if len(args) != 0:
    print "Error: not all arguments recognised\n"
    usage(); sys.exit()
if '-r' not in opts.keys():
    print "Error: -r argument missing\n"
    usage(); sys.exit()
if '-l' in opts.keys():
    ids_list = True
if '-q' in opts.keys():
    queries = True
if ids_list and queries:
    usage(); sys.exit()
if not ids_list and not queries:
    print "Error: provide taxonomy IDs with -l or -q\n"
    usage(); sys.exit()

Entrez.email = "j.vosseberg@uu.nl" # Identification required
ranks = opts['-r'].lower()
ranks = ranks.split(",")

if queries:
    tax_ids = opts['-q'].split(",")
    for tax_id in tax_ids:
        print "\t".join([tax_id] + get_rank(tax_id, ranks))
else:
    try:
        tax_ids_list = open(opts['-l'])
    except IOError:
        print opts['-l'], "not found"; sys.exit()
    for tax_id in tax_ids_list:
        tax_id = tax_id.rstrip()
        print "\t".join([tax_id] + get_rank(tax_id, ranks))
    tax_ids_list.close()
