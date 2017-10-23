from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse
import time
import logging

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Edited Martin Aryee's function to demultiplex to reduce
# the amount of memory required to run.
# Martin's original demultiplex function can be found
# here:
# https://github.com/aryeelab/umi/wiki
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__ = 'Allison MacLeay'

logging.basicConfig()
logger = logging.getLogger('demultiplex_nobuff')
logger.setLevel(logging.INFO)


def fq(file):
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')
    with fastq as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            yield [l1, l2, l3, l4]


def create_key(f1, f2):
    """Create a combined barcode key file from 2 files"""
    bc_dictA = {}
    bc_dictP = {}
    add_file_to_dict(f1, bc_dictA)
    add_file_to_dict(f2, bc_dictP)
    bcAP = [bc_dictA, bc_dictP]
    return bcAP


def add_file_to_dict(fname, d):
    """
    helper function - add a file to the dictionary
    :param fname: an array of strings with filenames to parse
    :param d: a dictionary (dict)
    :return:
    """
    HEADER = 1  # skip first line if equal to 1
    fh = open(fname, 'r')
    for line in fh:
        if HEADER == 1:
            HEADER = HEADER-1
            continue
        line = line.strip().split('\t')
        if len(line) < 2:
            continue
        [id, seq] = line
        d[seq[1:8]] = id
    fh.close()
    return

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# helper functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_sample_name(i1,i2,id_array):
    index1=i1[1][1:8]
    index2=i2[1][1:8]
    dictA=id_array[0]
    dictP=id_array[1]
    Aname=''
    Pname=''
    if index1 in dictP.keys():
        Pname=dictP[index1]
    else:
        Pname=index1
    if index2 in dictA.keys():
        Aname=dictA[index2]
    else:
        Aname=index2
    return (Aname + '_' + Pname)
    
def get_seq(i1, i2):
    seq1=i1[1]
    seq2=i2[1]
    return seq1[1:8] + seq2[1:8]

def demultiplex(read1, read2, index1, index2, p5_barcodes, p7_barcodes, out_dir, out_fname=None, min_reads=10000):

    # args = {'out_dir':'/PHShome/ma695/tmp', 'min_reads':10}
    # base = '/data/joung/sequencing_bcl/131007_M01326_0075_000000000-A6B33/Data/Intensities/BaseCalls'
    # args['read1'] = os.path.join(base, 'Undetermined_S0_L001_R1_001.fastq.gz')
    # args['read2'] = os.path.join(base, 'Undetermined_S0_L001_R2_001.fastq.gz')
    # args['index1'] = os.path.join(base, 'Undetermined_S0_L001_I1_001.fastq.gz')
    # args['index2'] = os.path.join(base, 'Undetermined_S0_L001_I2_001.fastq.gz')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    bar_dict = [p5_barcodes, p7_barcodes] if isinstance(p5_barcodes, dict) else create_key(p5_barcodes, p7_barcodes)

    fname = ''
    if out_fname:
        fname = args['out_fname'] + '_'

    outfiles_r1 = {}
    outfiles_r2 = {}
    outfiles_i1 = {}
    outfiles_i2 = {}

    total_count = 0
    count = {}

    # Create count dictionary first
    start = time.time()
    for i1, i2 in itertools.izip(fq(index1), fq(index2)):
        sample_id = get_seq(i1, i2)

        # Increment read count and create output buffers if this is a new sample barcode
        if not count.has_key(sample_id):
            count[sample_id] = 0
        count[sample_id] += 1
        total_count += 1
        if total_count % 5000000 == 0:
            logger.info("Processed %d counts in %.1f minutes." % (total_count, (time.time() - start) / 60))

    logger.info("Read count complete in %.1f minutes." % ((time.time() - start) / 60))

    total_count = 0
    for r1, r2, i1, i2 in itertools.izip(fq(read1), fq(read2), fq(index1), fq(index2)):
        # the original demultiplex stored sequences in a buffer to execute in 1N instead of 2N
        # this version minimizes the memory requirement by running in 2N
        total_count += 1
        if total_count % 1000000 == 0:
            logger.info("Processed %d reads in %.1f minutes." % (total_count, (time.time() - start) / 60))
        sample_id = get_sample_name(i1, i2, bar_dict)
        if count[get_seq(i1, i2)] < min_reads:
            # Write remaining buffered reads to a single fastq.
            # (These reads correspond to barcodes that were seen less than min_reads times)
            if 'undetermined_r1' not in vars():
                undetermined_r1 = open(os.path.join(out_dir, fname + 'undetermined.r1.fastq'), 'w')
            if 'undetermined_r2' not in vars():
                undetermined_r2 = open(os.path.join(out_dir, fname + 'undetermined.r2.fastq'), 'w')
            if 'undetermined_i1' not in vars():
                undetermined_i1 = open(os.path.join(out_dir, fname + 'undetermined.i1.fastq'), 'w')
            if 'undetermined_i2' not in vars():
                undetermined_i2 = open(os.path.join(out_dir, fname + 'undetermined.i2.fastq'), 'w')
            for line in r1:
                print(line, file=undetermined_r1, end="")
            for line in r2:
                print(line, file=undetermined_r2, end="")
            for line in i1:
                print(line, file=undetermined_i1, end="")
            for line in i2:
                print(line, file=undetermined_i2, end="")
        else:
            if sample_id not in outfiles_r1.keys():
                outname = fname + sample_id
                outfiles_r1[sample_id] = open(os.path.join(out_dir, '%s.r1.fastq' % outname), 'w')
                outfiles_r2[sample_id] = open(os.path.join(out_dir, '%s.r2.fastq' % outname), 'w')
                outfiles_i1[sample_id] = open(os.path.join(out_dir, '%s.i1.fastq' % outname), 'w')
                outfiles_i2[sample_id] = open(os.path.join(out_dir, '%s.i2.fastq' % outname), 'w')
            for line in r1:
                print(line, file=outfiles_r1[sample_id], end="")
            for line in r2:
                print(line, file=outfiles_r2[sample_id], end="")
            for line in i1:
                print(line, file=outfiles_i1[sample_id], end="")
            for line in i2:
                print(line, file=outfiles_i2[sample_id], end="")

    undetermined_r1.close()
    undetermined_r2.close()
    undetermined_i1.close()
    undetermined_i2.close()

    for sample_id in outfiles_r1.keys():
        outfiles_r1[sample_id].close()
    for sample_id in outfiles_r2.keys():
        outfiles_r2[sample_id].close()
    for sample_id in outfiles_i1.keys():
        outfiles_i1[sample_id].close()
    for sample_id in outfiles_i2.keys():
        outfiles_i2[sample_id].close()

    num_fastqs = len([v for k, v in count.iteritems() if v >= min_reads])
    logger.info('Wrote FASTQs for the %d sample barcodes out of %d with at least %d reads in %.1f minutes.' % (
    num_fastqs, len(count), min_reads, (time.time() - start) / 60))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#          MAIN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--read1', required=True)
    parser.add_argument('--read2', required=True)
    parser.add_argument('--index1', required=True)
    parser.add_argument('--index2', required=True)
    parser.add_argument('--min_reads', type=int, default=10000)
    parser.add_argument('--p5_barcodes')
    parser.add_argument('--p7_barcodes')
    parser.add_argument('--out_dir', default='.')
    parser.add_argument('--out_fname', default='')
    args = vars(parser.parse_args())
    swap = {}
    do_swap = 1
    fargs = ['read1', 'read2', 'index1', 'index2']
    for f in fargs:
        name = args[f]
        if name.find('_R1_') > 0:
            swap['read1'] = name
        elif name.find('_R2_') > 0:
            swap['read2'] = name
        elif name.find('_I1_') > 0:
            swap['index1'] = name
        elif name.find('_I2_') > 0:
            swap['index2'] = name
        else:
            do_swap = 0  # one or more files do not adhere to schema.  Can not confidently swap

    if (do_swap == 1) & (len(swap) == 4):
        # swapping files for names that are passed in in the wrong order
        # but follow the schema of containing _R1_
        for f in fargs:
            args[f] = swap[f]
    demultiplex(args['read1'], args['read2'], args['index1'], args['index2'], args['p5_barcodes'], args['p7_barcodes'], args['out_dir'], args['out_fname'],
                min_reads=args['min_reads'])
