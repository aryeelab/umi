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
Retrieve the fastq file (allowing gz formats) as 'yield' objects
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
Also will write a table with all the barcodes and total number of reads found in fastq file 01
"""
def fq_demultiplex(fastq_file01, fastq_file02, barcodes, min_reads, count_report, out_name01, out_name02, out_dir, prefix_table):
    # Comments about input values
    if barcodes == []:
        print('\nDemultiplexing without a barcode list, and with a minimum of %s reads' % min_reads, file=sys.stderr)
    else:
        print('\nDemultiplexing only samples with barcode from the input list', file=sys.stderr)
    if int(min_reads) <= 0 or int(count_report) <= 0:
        print('ERROR: numerical input values must be positive', file=sys.stderr)
        return()        

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # The buffers are used to avoid writing fastq files for sample barcodes with total reads less than min_reads
    outfiles_r01 = {}
    outfiles_r02 = {}

    buffer_r01 = {}
    buffer_r02 = {}

    count = {}
    bc_count = {}
    total_count = 0

    # Report counted samples every 'count_report'
    start = time.time()
    for item01, item02 in itertools.izip(fq_open(fastq_file01), fq_open(fastq_file02)):
        total_count += 1
        if total_count % float(count_report) == 0:
            print("Processed %d reads in %.1f minutes." % (total_count, (time.time()-start)/60), file=sys.stderr)

        sample_barcode = item01[0].split()[2]

        # Count sample barcode's read, independently if barcode is an input list or not
        if not bc_count.has_key(sample_barcode):
            bc_count[sample_barcode] = 0
        bc_count[sample_barcode] += 1

        # Create output buffer if this is a new sample barcode
        if not count.has_key(sample_barcode):
            count[sample_barcode] = 0
            buffer_r01[sample_barcode], buffer_r02[sample_barcode] = list(), list()

        # Increment read count according to whether or not a sample barcode list was given
        if len(barcodes) > 0:
             if sample_barcode in barcodes:
                count[sample_barcode] += 1
        else:
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

    # Write table with barcode and number of reads for every barcode
    with open(os.path.join(out_dir, '_'.join([prefix_table, 'Sample_Barcodes_Reads.txt'])), 'w') as out_bclist:
        print('barcode', 'total_reads', sep='\t', file=out_bclist)
        for key, value in bc_count.items():
            if value > 0:
                print(key, value, sep='\t', file=out_bclist)
 
    num_fastqs = len([v for k,v in count.iteritems() if v >= int(min_reads)] )
    print ('Wrote FASTQs for %d sample barcodes out of %d with at least %d read(s).' %(num_fastqs, len(count), int(min_reads)), file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Demultiplex fastq file w.r.t. sample barcode, which is the third field in read header')
    parser.add_argument('--fastq_file01', help='re-tagged fastq file with Reads 01', required=True)
    parser.add_argument('--fastq_file02', help='re-tagged fastq file with Reads 02', required=True)
    parser.add_argument('--barcodes', help='sample barcodes file without header', default=[])
    parser.add_argument('--min_reads', help='minimum number of reads required for a sample with given barcode to be stored', default=500)
    parser.add_argument('--count_report', help='counted samples will be reported every count_report', default=3000000)
    parser.add_argument('--out_name01', help='base name for output fastq associated to Read 1, including the ".fastq" file extension', required=True)
    parser.add_argument('--out_name02', help='base name for output fastq associated to Read 2, including the ".fastq" file extension', required=True)
    parser.add_argument('--out_dir', help='directory where to save output files', default='.')
    parser.add_argument('--prefix_table', help='prefix for table with all the barcodes and their number of reads associated', default='')
    args = parser.parse_args()

    barcodes = list()
    barcodes_list = args.barcodes
    if not isinstance(barcodes_list, list):
        with open(args.barcodes, 'r') as f:
            lines = f.readlines()
            for bc in lines:
                barcodes.append(bc.rstrip())

    fq_demultiplex(args.fastq_file01, args.fastq_file02, barcodes, args.min_reads, args.count_report, args.out_name01, args.out_name02, args.out_dir, args.prefix_table)
    

if __name__ == "__main__":
    main()
