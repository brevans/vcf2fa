vcf2fa
======

Given a vcf and a set of bam files, generate consensus sequences for each individual (taking into account areas of low/no coverage)

Required Software:
- python http://www.python.org/downloads/
- biopython http://biopython.org/wiki/Download
- samtools http://sourceforge.net/projects/samtools/files/
- bedtools https://github.com/arq5x/bedtools2


Prerequisites:
- reference fasta
- sorted, indexed bam files for each sample, mapped to your reference
- vcf file corresponding to bam files
 e.g.
```
samtools mpileup -uDf ref.fasta *.bam | bcftools view -vcg - > var.vcf
```
STEP 1:
```bash
#Generate a matrix of depth of coverage per sample per position:
./gen_bed_files.py reference.fa
bedtools multicov -bams path/to/bams/*.bam -bed reference_single_base.bed > multicov.bed
```
