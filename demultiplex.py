from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse

__author__ = 'Martin Aryee'


parser = argparse.ArgumentParser()
parser.add_argument('--read1', required=True)
parser.add_argument('--read2', required=True)
parser.add_argument('--index1', required=True)
parser.add_argument('--index2', required=True)
parser.add_argument('--min_reads', type=int, default=10000)
parser.add_argument('--out_dir', default='.')
args = vars(parser.parse_args())
out_dir = args['out_dir']

#args = {'out_dir':'/PHShome/ma695/tmp', 'min_reads':10}
#base = '/data/joung/sequencing_bcl/131007_M01326_0075_000000000-A6B33/Data/Intensities/BaseCalls'
#args['read1'] = os.path.join(base, 'Undetermined_S0_L001_R1_001.fastq.gz')
#args['read2'] = os.path.join(base, 'Undetermined_S0_L001_R2_001.fastq.gz')
#args['index1'] = os.path.join(base, 'Undetermined_S0_L001_I1_001.fastq.gz')
#args['index2'] = os.path.join(base, 'Undetermined_S0_L001_I2_001.fastq.gz')



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

def get_sample_barcode(i1, i2):
    seq1 = i1[1]
    seq2 = i2[1]
    return seq1[0:8] + seq2[0:8]

# Create molecular ID by concatenating molecular barcode and beginning of r1 read sequence
def get_umi(r1, r2, i1, i2):
    molecular_barcode = i2[1][8:16]
    return '%s_%s' % (molecular_barcode, r1[1][0:6])


if not os.path.exists(out_dir):
    os.makedirs(out_dir)

outfiles_r1 = {}
outfiles_r2 = {}
outfiles_i1 = {}
outfiles_i2 = {}

count = {}
buffer_r1 = {}
buffer_r2 = {}
buffer_i1 = {}
buffer_i2 = {}

#it = itertools.izip(fq(args['read1']), fq(args['read2']), fq(args['index1']), fq(args['index2']))
#for r1,r2,i1,i2 in itertools.islice(it, 0, 100):
for r1,r2,i1,i2 in itertools.izip(fq(args['read1']), fq(args['read2']), fq(args['index1']), fq(args['index2'])):
    sample_barcode = get_sample_barcode(i1, i2)
    # Create molecular ID by concatenating molecular barcode and beginning of r1 read sequence
    molecular_id = get_umi(r1, r2, i1, i2)
    # Add molecular id to read headers
    r1[0] = '%s %s\n' % (r1[0].rstrip(), molecular_id)
    r2[0] = '%s %s\n' % (r2[0].rstrip(), molecular_id)
    i1[0] = '%s %s\n' % (i1[0].rstrip(), molecular_id)
    i2[0] = '%s %s\n' % (i2[0].rstrip(), molecular_id)

    # Increment read count and create output buffers if this is a new sample barcode
    if not count.has_key(sample_barcode):
        count[sample_barcode] = 0
        buffer_r1[sample_barcode] = []
        buffer_r2[sample_barcode] = []
        buffer_i1[sample_barcode] = []
        buffer_i2[sample_barcode] = []
    count[sample_barcode] += 1

    # Write the reads to file or to sample-specific output buffers if we haven't
    # reached min_reads
    # The buffers are used to avoid writing fastqs for sample barcodes with very few reads.
    if count[sample_barcode] < args['min_reads']:
        buffer_r1[sample_barcode].append(r1)
        buffer_r2[sample_barcode].append(r2)
        buffer_i1[sample_barcode].append(i1)
        buffer_i2[sample_barcode].append(i2)
    elif count[sample_barcode] == args['min_reads']:
        outfiles_r1[sample_barcode] = open(os.path.join(out_dir, '%s.R1.fastq' % sample_barcode), 'w')
        outfiles_r2[sample_barcode] = open(os.path.join(out_dir, '%s.R2.fastq' % sample_barcode), 'w')
        outfiles_i1[sample_barcode] = open(os.path.join(out_dir, '%s.I1.fastq' % sample_barcode), 'w')
        outfiles_i2[sample_barcode] = open(os.path.join(out_dir, '%s.I2.fastq' % sample_barcode), 'w')
        # Spill the buffers to sample-specific fastqs
        for record in buffer_r1[sample_barcode] + r1:
            print (''.join(record), file=outfiles_r1[sample_barcode], end="")
        for record in buffer_r2[sample_barcode] + r2:
            print (''.join(record), file=outfiles_r2[sample_barcode], end="")
        for record in buffer_i1[sample_barcode] + i1:
            print (''.join(record), file=outfiles_i1[sample_barcode], end="")
        for record in buffer_i2[sample_barcode] + i2:
            print (''.join(record), file=outfiles_i2[sample_barcode], end="")
        del buffer_r1[sample_barcode]
        del buffer_r2[sample_barcode]
        del buffer_i1[sample_barcode]
        del buffer_i2[sample_barcode]
    else:
        for line in r1:
            print (line, file=outfiles_r1[sample_barcode], end="")
        for line in r2:
            print (line, file=outfiles_r2[sample_barcode], end="")
        for line in i1:
            print (line, file=outfiles_i1[sample_barcode], end="")
        for line in i2:
            print (line, file=outfiles_i2[sample_barcode], end="")

# Write remaining buffered reads to a single fastq.
# (These reads correspond to barcodes that were seen less than min_reads times)
undetermined_r1 = open(os.path.join(out_dir, 'undetermined.R1.fastq'), 'w')
undetermined_r2 = open(os.path.join(out_dir, 'undetermined.R2.fastq'), 'w')
undetermined_i1 = open(os.path.join(out_dir, 'undetermined.I1.fastq'), 'w')
undetermined_i2 = open(os.path.join(out_dir, 'undetermined.I2.fastq'), 'w')
for sample_barcode in buffer_r1.keys():
    for record in buffer_r1[sample_barcode]:
        print (''.join(record), file=undetermined_r1, end="")
    for record in buffer_r2[sample_barcode]:
        print (''.join(record), file=undetermined_r2, end="")
    for record in buffer_i1[sample_barcode]:
        print (''.join(record), file=undetermined_i1, end="")
    for record in buffer_i2[sample_barcode]:
        print (''.join(record), file=undetermined_i2, end="")

# Close files
for sample_barcode in outfiles_r1:
    outfiles_r1[sample_barcode].close()
    outfiles_r2[sample_barcode].close()
    outfiles_i1[sample_barcode].close()
    outfiles_i2[sample_barcode].close()
undetermined_r1.close()
undetermined_r2.close()
undetermined_i1.close()
undetermined_i2.close()

num_fastqs = len([v for k,v in count.iteritems() if v>=args['min_reads']])
print('Wrote FASTQs for the %d sample barcodes out of %d with at least %d reads.' % (num_fastqs, len(count), args['min_reads']))
