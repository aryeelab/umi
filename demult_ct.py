from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse
import time
import logging
from collections import Counter
import pandas as pd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Edited Martin Aryee's function to demultiplex to reduce
# the amount of memory required to run.
# Demultiplex based on sample ID and molecular barcode
# for increased granularity when processing cell free
# tumor dna.
# Martin's original demultiplex function can be found
# here:
# https://github.com/aryeelab/umi/wiki
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

__author__ = 'Allison MacLeay'

logging.basicConfig()
logger = logging.getLogger('demult_ct')
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
    bc_dict_a = {}
    bc_dict_p = {}
    add_file_to_dict(f1, bc_dict_a)
    add_file_to_dict(f2, bc_dict_p)
    bc_ap = [bc_dict_a, bc_dict_p]
    return bc_ap


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


def get_sample_name(i1, i2, id_array):
    """Get sample name from p5 and p7 barcode dictionary"""
    index1 = i1[1][1:8]
    index2 = i2[1][1:8]
    dict_a = id_array[0]
    dict_p = id_array[1]

    if index1 in dict_p.keys():
        pname = dict_p[index1]
    else:
        pname = index1
    if index2 in dict_a.keys():
        aname = dict_a[index2]
    else:
        aname = index2
    return '{}_{}'.format(aname, pname)


def get_seq(i1, i2):
    """Get sample id sequence """
    seq1 = i1[1]
    seq2 = i2[1]
    return seq1[1:8] + seq2[1:8]


def get_molecular_barcode(read1, read2, index1, index2, strategy='8B12H,,,'):
    """ Return molecular barcode and processed R1 and R2
        molecular barcode is first 8 bases, 12 handle
        :param strategy - B=barcode H=handle "R1_string,R2_string"
        :return mol_bc, read1, read2
    """
    barcode = ''
    r1_out = ''
    r2_out = ''
    for i, (read, pstr) in enumerate(zip([read1, read2, index1, index2], strategy.split(','))):
        pos = 0
        for group in re.findall('\d+\D+', pstr):
            btype = group[-1]
            num = int(group[:-1])
            if btype == 'B':  # Barcode
                barcode += read[1][pos:pos + num]
            elif btype == 'H':  # Handle - skip this
                pass
            pos += num
        if i == 0:
            r1_out = read[1][pos:]
        elif i == 1:
            r2_out = read[1][pos:]
    return barcode, r1_out, r2_out


def ct_to_dict(count):
    """Count data to dictionary """
    ctd = {'sample': [], 'molecular_barcode': [], 'count': []}
    for sample in count:
        for mol_bc in count[sample]:
            ctd['sample'].append(sample)
            ctd['molecular_barcode'].append(mol_bc)
            ctd['count'].append(count[sample][mol_bc])
    return ctd


def demultiplex(read1, read2, out_dir, index1=None, index2=None, p5_barcodes=None, p7_barcodes=None, out_fname=None,
                min_reads=10000, min_mol_bc=100, stats_out=None):
    """ Demultiplex based on sample name and molecular barcode """
    mode = 0  # 0 - demultiplex, 1 - molecular barcode only
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if None not in [p5_barcodes, p7_barcodes]:
        bar_dict = [p5_barcodes, p7_barcodes] if isinstance(p5_barcodes, dict) else create_key(p5_barcodes, p7_barcodes)
    else:
        logger.info('No p5 and p7 file provided.  Demultiplexing on molecular barcode only.')
        index1 = read1  # placeholder
        index2 = read2
        mode = 1

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
    for r1, r2, i1, i2 in itertools.izip(fq(read1), fq(read2), fq(index1), fq(index2)):
        sample_id = get_seq(i1, i2) if mode == 0 else os.path.basename(read1).split('.')[0]
        mol_bc, _, _ = get_molecular_barcode(r1, r2, i1, i2)

        # Increment read count and create output buffers if this is a new sample barcode
        if not count.has_key(sample_id):
            count[sample_id] = Counter()
        count[sample_id][mol_bc] += 1

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
        sample_index = get_seq(i1, i2) if mode == 0 else os.path.basename(read1).split('.')[0]
        sample_id = get_sample_name(i1, i2, bar_dict) if mode == 0 else sample_index
        mol_bc, r1_seq, r2_seq = get_molecular_barcode(r1, r2, i1, i2)
        if sum(count[sample_index].values()) < min_reads and mode == 0:
            # Write remaining buffered reads to a single fastq.
            # (These reads correspond to barcodes that were seen less than min_reads times)
            sample_id = 'undetermined'
            mol_bc = 'molbc'
        elif count[sample_index][mol_bc] < min_mol_bc:
            mol_bc = 'molbc'
        sample_id += '_{}'.format(mol_bc)

        if sample_id not in outfiles_r1.keys():
            outname = fname + sample_id
            outfiles_r1[sample_id] = open(os.path.join(out_dir, '%s.r1.fastq' % outname), 'w')
            outfiles_r2[sample_id] = open(os.path.join(out_dir, '%s.r2.fastq' % outname), 'w')
            if mode == 0:
                outfiles_i1[sample_id] = open(os.path.join(out_dir, '%s.i1.fastq' % outname), 'w')
                outfiles_i2[sample_id] = open(os.path.join(out_dir, '%s.i2.fastq' % outname), 'w')
        for line in r1:
            print(line, file=outfiles_r1[sample_id], end="")
        for line in r2:
            print(line, file=outfiles_r2[sample_id], end="")
        if mode == 0:
            for line in i1:
                print(line, file=outfiles_i1[sample_id], end="")
            for line in i2:
                print(line, file=outfiles_i2[sample_id], end="")

    for sample_id in outfiles_r1.keys():
        outfiles_r1[sample_id].close()
    for sample_id in outfiles_r2.keys():
        outfiles_r2[sample_id].close()
    if mode == 0:
        for sample_id in outfiles_i1.keys():
            outfiles_i1[sample_id].close()
        for sample_id in outfiles_i2.keys():
            outfiles_i2[sample_id].close()

    num_fastqs = len([v for k, v in count.iteritems() if sum(v.values()) >= min_reads])
    logger.info('Wrote FASTQs for the %d sample barcodes out of %d with at least %d reads in %.1f minutes.' % (
    num_fastqs, len(count), min_reads, (time.time() - start) / 60))
    if stats_out:
        df = pd.DataFrame(ct_to_dict(count))
        made_cut = pd.DataFrame(df.groupby('sample')['count'].sum() > min_reads)
        mdf = pd.merge(df, made_cut, left_on='sample', right_index=True, suffixes=['', '_sample_pass'])
        mdf['molecular_bc_passed'] = mdf['count'] > min_mol_bc
        mdf['written'] = mdf['count_sample_pass'] & mdf['molecular_bc_passed']
        mdf.to_csv(stats_out, sep='\t', index=False)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#          MAIN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--read1', required=True)
    parser.add_argument('--read2', required=True)
    parser.add_argument('--index1')  # if not provided split by molecular barcode and not sample id
    parser.add_argument('--index2')
    parser.add_argument('--min_reads', type=int, default=10000)
    parser.add_argument('--min_barcodes', type=int, default=100)
    parser.add_argument('--p5_barcodes')
    parser.add_argument('--p7_barcodes')
    parser.add_argument('--out_dir', default='.')
    parser.add_argument('--out_fname', default='')
    parser.add_argument('--stats_out', default=None)
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
    demultiplex(args['read1'], args['read2'], args['out_dir'], args['index1'], args['index2'], args['p5_barcodes'],
                args['p7_barcodes'], args['out_fname'], min_reads=args['min_reads'], min_mol_bc=args['min_barcodes'],
                stats_out=None)
