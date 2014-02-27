#!/usr/bin/env python
from os import path
from os import makedirs
from collections import defaultdict as dd
import re
#from Bio import SeqIO

ambi = {'AC':'M', 'AG':'R', 'AT':'W', 'CG':'S', 'CT':'Y', 'GT':'K',
        'CA':'M', 'GA':'R', 'TA':'W', 'GC':'S', 'TC':'Y', 'TG':'K',
        'AA':'A', 'CC':'C', 'GG':'G', 'TT':'T'}

def m_dir(d):
    try:
        makedirs(d)
    except OSError:
        pass

class rec(object):
    def __init__ (self):
        '''
        container for previous, current, next objects
        '''
        self.up=0
        self.curr=0
        self.down=0

class Vcf(object):
    '''
    fVcf: POSITION SORTED vcf file name
    vcf: a generator that walks over a vcf file

    sample VCF line:
    CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT sam1...samn
    '''
    def __init__(self, fVcf):
        self.fVcf = fVcf
        self.header_lines=[]
        self.vcf_handle=open(fVcf)
        self.rec = rec()

        l=''
        #initialize samples, strip off vcf header
        #stop at first line with real data
        while True:
            l = self.vcf_handle.readline()

            if l.startswith("#CHROM"):
                a=l.lstrip('#').rstrip().split('\t')
                self.cols = a[:9]
                self.samples = []
                for sam in a[9:]:
                    sam_name = re.sub('\.bam$|\.sorted\.bam$', '', path.basename(sam))
                    self.samples.append(sam_name)

            elif l.startswith("##"):
                self.header_lines.append(l.rstrip())

            else:
                self.rec.up = 0
                self.rec.curr = 0
                self.rec.down = self.parse_vcf_line(l)
                self.format = self.rec.down['FORMAT'].split(':')
                break
    
    def __iter__(self):
        return self

    def next(self):
        '''
        '''
        if self.rec.curr != None:
            self.rec.up = self.rec.curr
            self.rec.curr = self.rec.down
            self.rec.down = self.parse_vcf_line(self.vcf_handle.readline())
            return self.rec.curr
        else:
            raise StopIteration()

    def parse_vcf_line(self, line):
        tmp = {}
        vals = line.split()
        if vals[7].startswith('INDEL'):
            vals = self.vcf_handle.readline().split()
        if len(vals) != len(self.cols)+len(self.samples):
            return None
        for i,v in zip(self.cols,vals[:9]):
            tmp[i]=v
        for i,v in zip(self.samples, vals[9:]):
            try:
                tmp[i] = dict(zip(self.format, v.split(':')))
            except AttributeError:
                tmp[i] = dict(zip(tmp['FORMAT'].split(':'), v.split(':')))
        return tmp

def vcf_gt_to_dna(snp, ind):
    ref = snp['REF']
    alt = snp['ALT']
    gt = snp[ind]['GT']
    poss_gts = [ref.upper()]+[x.upper() for x in alt.split(',')] 
    trans = dict(zip(range(len(poss_gts)), poss_gts))
    sam_dna = ''.join([trans[x] for x in gt.split('/')])
    return sam_dna

def make_fasta(mincov, mcb_fn, vcf_fn):
    vcf = Vcf(vcf_fn)
    mcb = open(mcb_fn)
    m_dir('consensus')
    seqs = dd(lambda: [])

    curr_chrom = ''
    #get a snp
    try:
        snp = vcf.next()
    except StopIteration:
        snp = None

    for l in mcb:
        
        #get info about current position
        tmp = l.split()
        chrom = tmp[0]
        if chrom != curr_chrom:
            if len(seqs) != 0:
                out = open(path.join('consensus', '{}_consensus.fa'.format(curr_chrom)), 'w')
                for s in seqs.keys():
                    #write multi fasta
                    out.write('>{}\n{}\n'.format(s, ''.join(seqs[s])))
                out.close()
            seqs = dd(lambda: [])
            curr_chrom = chrom
        (pos0, pos1) = [int(x) for x in tmp[1:3]]
        ref_bp = tmp[3]
        cov = [int(x) for x in tmp[4:]]

        if snp == None or pos1 < snp['POS'] or pos1 > snp['POS']:
            for s,c in zip(vcf.samples,cov):
                if c >= mincov:
                    seqs[s].append(ref_bp)
                else:
                    seqs[s].append('N')
        elif pos1 == snp['POS']:
            for s,c in zip(vcf.samples,cov):
                if c >= mincov:
                    seqs[s].append(vcf_gt_to_dna(snp, s))
                else:
                    seqs[s].append('N')
            try:
                snp = vcf.next()
            except StopIteration:
                snp = None
    
    #write last chrom's worth of data
    out = open(path.join('consensus', '{}_consensus.fa'.format(curr_chrom)), 'w')
    for s in seqs.keys():
        #write multi fasta
        out.write('>{}\n{}\n'.format(s, ''.join(seqs[s])))
    out.close()
    

if __name__ == '__main__':
    make_fasta(7, 'multicov.bed', 'var.vcf')