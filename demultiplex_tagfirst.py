"""
demultiplex_tagfirst.py
"""

from __future__ import print_function

___author___ = "Martin Aryee"

import gzip
import re
import sys
import os
import argparse
import time
import itertools

"""
Retrieve the fastq file (allowing gz formats) as 'yield-generators' objects
"""
def fq_open(fastq_file):
    if re.search(".gz$", fastq_file):
        fastq = gzip.open(fastq_file, "rb")
    else:
        fastq = open(fastq_file, "r")

    with fastq as f:
        while True:
            line1 = f.readline()
            if not line1:
                break
            line2 = f.readline()
            line3 = f.readline()
            line4 = f.readline()
            yield [line1, line2, line3, line4]


"""
Demultiplex for re-tagged fastq files allocating by Sample Barcode
"""
def fq_demultiplex(fastq_file01, fastq_file02, barcodes, min_reads, count_report, out_name01, out_name02, out_dir):
    # Comments about input values
    if int(min_reads) <= 0 or int(count_report) <= 0:
        print('ERROR: numerical input values must be positive', file=sys.stderr)
        return ()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # The buffers are used to avoid writing fastq files for sample barcodes with total reads less than min_reads
    outfiles_r01, outfiles_r02 = {}, {}

    buffer_r01, buffer_r02 = {}, {}

    count = {}
    total_count = 0

    # Report counted samples every 'count_report'
    start = time.time()
    for item01, item02 in itertools.izip(fq_open(fastq_file01), fq_open(fastq_file02)):
        total_count += 1
        if total_count % float(count_report) == 0:
            print("Processed %d reads in %.1f minutes." % (total_count, (time.time() - start) / 60), file=sys.stderr)

        sample_barcode = item01[0].split()[2]

        # Create output buffer if this is a new sample barcode
        if sample_barcode not in count:
            count[sample_barcode] = 0
            buffer_r01[sample_barcode], buffer_r02[sample_barcode] = list(), list()

        # Increment read count according to whether or not the sample barcode is in the input barcode list
        if sample_barcode in barcodes:
            count[sample_barcode] += 1

        # Write reads to file or to specific buffer if we have not reached min_reads
        if count[sample_barcode] < int(min_reads):
            buffer_r01[sample_barcode].append(item01)
            buffer_r02[sample_barcode].append(item02)
        elif count[sample_barcode] == int(min_reads):
            outfiles_r01[sample_barcode] = open(os.path.join(out_dir, '_'.join([sample_barcode, out_name01])), 'w')
            for record01 in buffer_r01[sample_barcode] + item01:
                outfiles_r01[sample_barcode].write(''.join(record01))
            del buffer_r01[sample_barcode]
            outfiles_r02[sample_barcode] = open(os.path.join(out_dir, '_'.join([sample_barcode, out_name02])), 'w')
            for record02 in buffer_r02[sample_barcode] + item02:
                outfiles_r02[sample_barcode].write(''.join(record02))
            del buffer_r02[sample_barcode]
        else:
            for line in item01:
                print(line, file=outfiles_r01[sample_barcode], end="")
            for line in item02:
                print(line, file=outfiles_r02[sample_barcode], end="")

    # Write non_allocated buffered reads to a single fastq file
    undetermined_r01 = open(os.path.join(out_dir, '_'.join(['NON_allocated', out_name01])), 'w')
    undetermined_r02 = open(os.path.join(out_dir, '_'.join(['NON_allocated', out_name02])), 'w')

    for bc in buffer_r01.keys():
        for record in buffer_r01[bc]:
            undetermined_r01.write(''.join(record))
        for record in buffer_r02[bc]:
            undetermined_r02.write(''.join(record))
    for bc in outfiles_r01.keys():
        outfiles_r01[bc].close()
        outfiles_r02[bc].close()
    undetermined_r01.close()
    undetermined_r02.close()


def main():
    parser = argparse.ArgumentParser(description='Demultiplex fastq file w.r.t. sample barcode, which is the third field in read header')
    parser.add_argument('--fastq_file01', help='re-tagged fastq file with Reads 01', required=True)
    parser.add_argument('--fastq_file02', help='re-tagged fastq file with Reads 02', required=True)
    parser.add_argument('--barcodes', help='sample barcodes file without header', required=True)
    parser.add_argument('--min_reads', help='minimum number of reads required for a sample with given barcode to be stored', default=500)
    parser.add_argument('--count_report', help='counted samples will be reported every count_report', default=3000000)
    parser.add_argument('--out_name01', help='base name for output fastq associated to Read 1, including the ".fastq" file extension', required=True)
    parser.add_argument('--out_name02', help='base name for output fastq associated to Read 2, including the ".fastq" file extension', required=True)
    parser.add_argument('--out_dir', help='directory where to save output files', default='.')
    args = parser.parse_args()

    barcodes = list()
    barcodes_list = args.barcodes
    if not isinstance(barcodes_list, list):
        with open(args.barcodes, 'r') as f:
            for bc in f:
                barcodes.append(bc.rstrip())

    fq_demultiplex(args.fastq_file01, args.fastq_file02, barcodes, args.min_reads, args.count_report, args.out_name01, args.out_name02, args.out_dir)


if __name__ == "__main__":
    main()
