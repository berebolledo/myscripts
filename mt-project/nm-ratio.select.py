#! /usr/bin/python
import sys, pysam


def get_nm(read):
    tags=dict(read.tags)
    return float(tags['NM'])


sample=sys.argv[1]
sam=pysam.Samfile(sample,'rb')
out=pysam.Samfile("nm-ratio."+sample,'wb',template=sam)


for read in sam:
    try:
        read1=read
        read2=sam.next()
        if get_nm(read1)<=(0.02*read1.rlen) and get_nm(read2)<=(0.02*read2.rlen):
            out.write(read1)
            out.write(read2)
        else:
            pass

    except StopIteration:
            sam.close()
            out.close()
