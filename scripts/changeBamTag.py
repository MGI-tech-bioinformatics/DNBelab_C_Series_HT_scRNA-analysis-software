import argparse,pysam
parser = argparse.ArgumentParser()
parser.add_argument('-i','--infile', metavar='FILE',default='anno_decon_sorted.bam')
parser.add_argument('-o','--outfile', metavar='FILE',default='anno_decon_CB.bam')
args = parser.parse_args()



import sys
import pysam
samfile = pysam.AlignmentFile(args.infile, "rb")
outsam = pysam.AlignmentFile(args.outfile, "wb",header=dict(samfile.header))

for i in samfile:
    try:
        if i.has_tag('UB'):
            RCB = i.get_tag('CB')
            RDB = i.get_tag('DB')
            i.set_tag('DB',RCB)
            i.set_tag('CB',RDB)
            outsam.write(i)
    except KeyError:
        continue

outsam.close()