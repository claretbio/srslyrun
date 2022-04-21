#!/usr/bin/env python3

import pysam
import argparse

def write_mapped_mol_len(bam, outfile):
    bamfile = pysam.AlignmentFile(bam, 'rb')
    with open(outfile, "w") as f:
        for b in bamfile:
            insertsize = b.template_length
            if insertsize != 0:
                f.write("{}\n".format(insertsize))

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("bamfile")
    ap.add_argument("outfile")
    a = ap.parse_args()
    write_mapped_mol_len(a.bamfile, a.outfile)
