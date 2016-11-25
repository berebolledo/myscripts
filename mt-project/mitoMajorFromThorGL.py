#! /usr/bin/python

# calculate the major allele, and extract the frequency of the non-major allele
# from Thor's angsd output

import sys


def get_major(frequencies):
    alleles={0:'A',1:'C',2:'G',3:'T'}
    if sum(frequencies)==0:
        major_allele="N"
    else:
        idx=frequencies.index(max(frequencies))
        major_allele=alleles[idx]
    return major_allele




fasta=open(sys.argv[1]+'.major.fa','w+')
seq_list=['N']*16569

with open(sys.argv[1]) as data:
    for d in data:
        f=[float(n) for n in d.strip().split('\t')[2:6]]
        non_major=d.strip().split('\t')[6]
        pos=d.strip().split('\t')[1]
        major_allele=get_major(f)
        seq_list[int(pos)-1]=major_allele

seq=''.join(seq_list)
fasta.write('>chrM\n')

for i in range(0,len(seq),100):
   fasta.write(seq[i:i+100]+'\n')

fasta.close()
