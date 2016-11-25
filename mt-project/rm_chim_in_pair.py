#! /usr/bin/python

import sys,pysam


def check_chim(read):
    tags = dict(read.tags)
    if 'SA' in tags.keys():
        sa = tags['SA'].split(';')[:-1]
        if len(sa)==1:
            chrom,pos,strand,cigar,mapq,nm=sa[0].split(',')
            if chrom == 'chrM' and ( int(pos) >=16000 or int(pos) <=600 ):
                return read
            else:
                pass
        else:
            pass
    else:
        return read

sample=sys.argv[1]
readLength=int(float(sys.argv[2]))

sam=pysam.Samfile(sample,'rb')
out=pysam.Samfile("dechim.rlen."+sample,'wb',template=sam)
for read in sam:
    try:
        read1=read
        read2=sam.next()
        if check_chim(read1) and check_chim(read2) and\
           read1.rlen>=readLength and read2.rlen>=readLength:
            out.write(read1)
            out.write(read2)
        else:
            pass
    except StopIteration:
            sam.close()
            out.close()

