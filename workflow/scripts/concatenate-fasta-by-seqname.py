#!/bin/env python

import sys

seqs = {}
with open(sys.argv[1], "r") as fh:
    curr = ""
    for line in fh:
        if line.startswith(">"):
            curr = line
            if curr not in seqs:
                seqs[curr] = ""
        else:
            seqs[curr] += line.strip()

for ident, seq in seqs.items():
    print("{}{}".format(ident, seq))
