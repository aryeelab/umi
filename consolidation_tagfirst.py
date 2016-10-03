"""
consolidation_tagfirst.py
"""

from __future__ import print_function


___author___ = "Martin Aryee"


import HTSeq
import subprocess
import sys
import os
import argparse


"""
Sort input fastq file by molecular id
"""
def fq_sort_fastq(fastq_file):
    if not os.path.exists(''.join([os.path.dirname(fastq_file), '/sorted/'])):
        os.makedirs(''.join([os.path.dirname(fastq_file), '/sorted/']))
    command_line = 'cat %s | paste - - - - | sort -k4,4 -k1,1 | tr "\t" "\n" > %s' % (fastq_file, ''.join([os.path.dirname(fastq_file), '/sorted/', 'sorted_', os.path.basename(fastq_file)]))
    subprocess.check_call(command_line, shell=True, env=os.environ.copy())

"""
Distribute reads according to sample_id and molecular_id
"""
def fq_read_bins(fastq_file_sorted):
    infile = HTSeq.FastqReader(fastq_file_sorted)
    read_num = 0
    bin_reads = []
    cur_molecular_id = ''
    for read in infile:
        read_num += 1
        read_name, sample_id, sample_barcode, molecular_id = read.name.split(' ')
        if molecular_id == cur_molecular_id:
            bin_reads.append(read)
        else:
            if cur_molecular_id != '':
                yield cur_molecular_id, cur_sample_id, bin_reads
            cur_molecular_id = molecular_id
            cur_sample_id = sample_id
            bin_reads = [read]
    yield cur_molecular_id, cur_sample_id, bin_reads #the last bin

"""
Get most common base for adequate frequency
"""
def fq_consolidate_position(bases, quals, min_qual, min_freq):
    num, qual = {}, {}

    num['A'] = num['C'] = num['G'] = num['T'] = num['N'] = 0
    qual['A'] = qual['C'] = qual['G'] = qual['T'] = qual['N'] = 0
    for bb, qq in zip(bases, quals):
        if qq > float(min_qual):
            num[bb] += 1
        if qq > qual[bb]:
            qual[bb] = qq
    most_common_base = max(num.iterkeys(), key=(lambda key: num[key]))
    freq = float(num[most_common_base]) / len(bases)
    if freq > float(min_freq):
        return True, most_common_base, qual[most_common_base]
    else:
        return False, 'N', 0


"""
Consolidate reads by molecular id
"""
def fq_consolidate(fastq_file_sorted, min_qual, min_freq, out_dir, out_name):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    outfile = open(os.path.join(out_dir, out_name), 'w')
    bins = fq_read_bins(fastq_file_sorted)

    num_input_reads = 0
    num_consolidated_reads = 0
    num_successes = 0 #bases with successful consolidation
    num_bases = 0
    for cur_molecular_id, cur_sample_id, reads in bins:
        num_input_reads += len(reads)
        num_consolidated_reads += 1

        # Get all the bases and quals in the read
        read_bases = zip(*[list(read.seq) for read in reads])
        read_quals = zip(*[list(read.qual) for read in reads])

        # Iterate position by position
        consolidation_sucess, cons_seq, cons_qual = zip(*[fq_consolidate_position(bases, quals, min_qual, min_freq) for bases, quals in zip(read_bases, read_quals)])

        # Count consolidation successes and failures
        num_successes += sum(consolidation_sucess)
        num_bases += len(consolidation_sucess)

        # Write consolidated FASTQ read
        outfile.write('@%s_%d %s\n' % (cur_molecular_id, len(reads), cur_sample_id)) #Header: Molecular id, number of reads, 2nd incoming header field (includes sample id)
        outfile.write(''.join(cons_seq) + '\n')
        outfile.write('+ \n')
        outfile.write(''.join([chr(q + 33) for q in cons_qual]) + '\n')

    print("Read %d input reads" % num_input_reads, file=sys.stderr)
    print("Wrote %d consolidated reads" % num_consolidated_reads, file=sys.stderr)
    print("Successfully consolidated %d bases out of %d (%.2f%%)" % (num_successes, num_bases, 100 * float(num_successes) / num_bases), file=sys.stderr)
    outfile.close()


def main():
    parser = argparse.ArgumentParser(description='Consolidate samples from a fastq file according to their molecular id. Run each paired file SEPARATELY')

    parser.add_argument('--fastq_file', help='fastq file', required=True)
    parser.add_argument('--min_qual', help='quality threshold', default=15)
    parser.add_argument('--min_freq', help='frequency threshold', default=0.8)
    parser.add_argument('--sorted', help='indicates if the input fastq file is already sorted', dest='sorted', action='store_true', default=False)
    parser.add_argument('--out_dir', help='directory where to save output files', default='.')
    parser.add_argument('--out_name', help='output fa=ile name, including extension ', required=True)

    args = parser.parse_args()

    if args.sorted:
        fastq_file_sorted = args.fastq_file
    else:    
        fq_sort_fastq(args.fastq_file)
        fastq_file_sorted = os.path.join(''.join([os.path.dirname(args.fastq_file), '/sorted/']), ''.join(['sorted_', os.path.basename(args.fastq_file)]))
    
    fq_consolidate(fastq_file_sorted, args.min_qual, args.min_freq, args.out_dir, args.out_name)
    

if __name__ == "__main__":
    main()
