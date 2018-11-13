#!/usr/bin/env python


from Bio import SeqIO
from optparse import OptionParser
# import numpy as np
import sys
import networkx as nx
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-i", "--hits", dest="hitspath", default="None",
                  help="[Required] Location of the HMMER hits fullseq fasta file")
parser.add_option("-t", "--taxadb", dest="taxapath", default="None",
                  help="[Required] Location of the file containing the taxonomy database from uniprot generated by extract_taxonomy.py")


# Parse options into variables
(options, args) = parser.parse_args()

hitspath = options.hitspath
taxapath = options.taxapath
if hitspath is None or taxapath is None:
    print "Invalid options"
    sys.exit(1)


Taxonomy_Tree = nx.DiGraph()

taxaDB = {}
f = open(taxapath)
lines = f.readlines()
for line in lines:
    ID, taxa = line.strip().split('\t')
    taxa = taxa.split(',')
    taxaDB[ID] = taxa


handle = SeqIO.parse(hitspath, 'fasta')
hits_sorted = {}
for Seq in handle:
    # ID = Seq.name.strip().split('_')[1]
    ID = Seq.name.strip().split('|')[0].split('_')[1]
    evalue = float(Seq.name.strip().split('|')[2])
    Seq.description = Seq.name.strip().split('|')[0]
    edges = []
    for t in range(len(taxaDB[ID])):
        Taxonomy_Tree.add_node(taxaDB[ID][t])
        if t > 0:
            edges.append((taxaDB[ID][t - 1], taxaDB[ID][t]))
    Taxonomy_Tree.add_node(ID)
    edges.append((taxaDB[ID][t], ID))
    Taxonomy_Tree.add_edges_from(edges)


nx.draw_networkx(Taxonomy_Tree, with_labels=True)
plt.show()
