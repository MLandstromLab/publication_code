import pysam
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Script for calculating CH score')
parser.add_argument('-v','--vcf', dest='vcf', help='VCF file')
parser.add_argument('-b','--bam', dest='bam', help='Bam file')
parser.add_argument('-o','--output', dest='output', help='Output file')
args = parser.parse_args()

# Read in the bam-file
samfile = pysam.AlignmentFile(args.bam, "rb")

def calculateCHscore(samfile_obj, chrom, pos, alt):

    fragment_lengths_alt = []
    for pileupcolumn in samfile_obj.pileup(chrom, pos, pos + 1):
        for pileupread in pileupcolumn.pileups:

            # Fragment length 
            tlen = pileupread.alignment.template_length

            # The base 
            base = pileupread.alignment.query_sequence[pileupread.query_position]

            if base == alt:
                fragment_lengths_alt.append(tlen)
            else:
                pass

    # Calculate the score 
    # 127–141 bp and 272–292 bp => tumor
    # 173–191 bp and 346–361 => CH
    tumor_bin_count = 0 
    ch_bin_count = 0
    for frag in fragment_lengths_alt:
        if ((127 =< fragment) & (fragment <= 141)):
            tumor_bin_count+=1 
        elif ..

    



