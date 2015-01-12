from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse

__author__ = 'Martin Aryee'

# python ~/work/projects/umi/consolidate.py --read1_in bcl/Undetermined_S0_L001_R1_001.fastq.gz --read2_in bcl/Undetermined_S0_L001_R2_001.fastq.gz --read1_out umi1.fastq.gz --read2_out umi2.fastq.gz --index1 bcl/Undetermined_S0_L001_I1_001.fastq.gz --index2 bcl/Undetermined_S0_L001_I2_001.fastq.gz


parser = argparse.ArgumentParser()
parser.add_argument('--read1_in', required=True)
parser.add_argument('--read2_in', required=True)
parser.add_argument('--read1_out', required=True)
parser.add_argument('--read2_out', required=True)
parser.add_argument('--index1', required=True)
parser.add_argument('--index2', required=True)
parser.add_argument('--out_dir', default='.')
args = vars(parser.parse_args())
out_dir = args['out_dir']

# args = {'out_dir':'/PHShome/ma695/tmp', 'min_reads':10}
# base = '/data/joung/sequencing_bcl/131007_M01326_0075_000000000-A6B33/Data/Intensities/BaseCalls'
# args['read1_in'] = os.path.join(base, 'Undetermined_S0_L001_R1_001.fastq.gz')
# args['read2_in'] = os.path.join(base, 'Undetermined_S0_L001_R2_001.fastq.gz')
# args['read1_out'] = os.path.join(out_dir, 'Undetermined_S0_L001_R1_001.umi.fastq.gz')
# args['read2_out'] = os.path.join(out_dir, 'Undetermined_S0_L001_R2_001.umi.fastq.gz')
# args['index1'] = os.path.join(base, 'Undetermined_S0_L001_I1_001.fastq.gz')
# args['index2'] = os.path.join(base, 'Undetermined_S0_L001_I2_001.fastq.gz')


def fq(file):
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')
    with fastq as f:
        while True:
            l1 = f.readline()
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            yield [l1, l2, l3, l4]


# Create molecular ID by concatenating molecular barcode and beginning of r1 read sequence
def get_umi(r1, r2, i1, i2):
    molecular_barcode = i2[1][8:16]
    return '%s_%s' % (molecular_barcode, r1[1][0:6])


if not os.path.exists(out_dir):
    os.makedirs(out_dir)

r1_out = gzip.open(args['read1_out'], 'wb')
r2_out = gzip.open(args['read2_out'], 'wb')

#for r1,r2,i1,i2 in itertools.izip(fq(args['read1']), fq(args['read2']), fq(args['index1']), fq(args['index2'])):
it = itertools.izip(fq(args['read1_in']), fq(args['read2_in']), fq(args['index1']), fq(args['index2']))
for r1,r2,i1,i2 in itertools.islice(it, 0, 100):
    # Create molecular ID by concatenating molecular barcode and beginning of r1 read sequence
    molecular_id = get_umi(r1, r2, i1, i2)
    # Add molecular id to read headers
    r1[0] = '%s %s\n' % (r1[0].rstrip(), molecular_id)
    r2[0] = '%s %s\n' % (r2[0].rstrip(), molecular_id)
    for line in r1:
        r1_out.write(line)
    for line in r2:
        r2_out.write(line)

# Close files
r1_out.close()
r2_out.close()


def sortFastq(fastq_toSort_dir):

    print("Sorting fastqs based on molecular barcode")

    for file in os.listdir(fastq_toSort_dir):
        if fnmatch.fnmatch(file, '*R1.fastq') or fnmatch.fnmatch(file, '*R2.fastq'):
            file_sorted = re.sub('.fastq$', '_sorted.fastq', file)
            cmd = 'cat ' + os.path.join(fastq_toSort_dir, file) + ' | paste - - - - | sort -k3,3 -t " " | tr "\t" "\n" >' + os.path.join(fastq_toSort_dir, file_sorted)
            subprocess.check_call(cmd, shell=True)
            os.remove(os.path.join(fastq_toSort_dir, file))
