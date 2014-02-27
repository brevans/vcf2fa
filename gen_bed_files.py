#!/usr/bin/env python
import sys
import re
from os import path
from Bio import SeqIO

ref_fn = sys.argv[1]
ref = SeqIO.index(ref_fn, "fasta")
sbed_fn = re.sub('\.fa$|\.fasta$', '_single_base.bed', path.basename(ref_fn))
genome_fn = re.sub('\.fa$|\.fasta$', '.genome', path.basename(ref_fn))
bed_fn = re.sub('\.fa$|\.fasta$', '.bed', path.basename(ref_fn))

with open(sbed_fn, 'w') as sbed_fh:
    with open(genome_fn, 'w') as genome_fh:
        with open(bed_fn, 'w') as bed_fh:
            for ch in ref.keys():
                genome_fh.write('{}\t{}\n'.format(ch, len(ref[ch])))
                bed_fh.write('{}\t0\t{}\n'.format(ch, len(ref[ch])))
                for pos in range(len(ref[ch])):
                    sbed_fh.write('{}\t{}\t{}\t{}\n'.format(
                                  ch, pos, pos+1, ref[ch][pos]))
