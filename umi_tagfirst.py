"""
umi_tagfirst.py
"""

from __future__ import print_function

___author___ = 'Martin Aryee'


import gzip
import re
import sys
import os
import argparse
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
Retrieve barcode and molecular_id of a read as two strings separated by an empty space if:
a) the EXACT restriction site is in the sequence;
b) the UMI tag and Barcode have FULL LENGTH.
Create molecular_id by concatenating UMI and beginning of r1 and r2 read sequences
"""
def fq_get_bcmolid(read01, read02, umi_len, restriction_enzyme_index, rs_len, molid_len):
    umi = read01[1][:int(umi_len)]
    barcode = read01[1][int(umi_len):int(restriction_enzyme_index)]
    id_pos = int(restriction_enzyme_index) + int(rs_len)
    molecular_id = '_'.join([ umi, read01[1][id_pos:(id_pos + int(molid_len))], read02[1][:int(molid_len)] ])
    out = '%s %s' % (barcode, molecular_id)
    return out


"""
Re-tag first line in the fastq files: seq_id line_id barcode umi:seq
Scores of barcode and umi are removed from line 4 in r1
"""
def fq_umitag(fastq_file01, fastq_file02, umi_len, barcode_len, restriction_enzyme, molid_len, out_name01, out_name02, out_dir):
    if int(barcode_len)<=0 or int(molid_len)<=0 or int(umi_len)<=0:
        print('ERROR: numerical input values must be positive', file=sys.stderr)
        return()        
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    out_path_name01 = os.path.join(out_dir, out_name01)
    out_path_name02 = os.path.join(out_dir, out_name02)
    non_out_path_name01 = os.path.join(out_dir, '_'.join(['NON_stored', out_name01]))
    non_out_path_name02 = os.path.join(out_dir, '_'.join(['NON_stored', out_name02]))

    # Files where we will write adequate elements
    fw = open(out_path_name01, "w")
    bw = open(out_path_name02, "w")

    # Files where we will write non-adequate elements
    non_fw = open(non_out_path_name01, "w")
    non_bw = open(non_out_path_name02, "w")

    rs_len = len(restriction_enzyme)

    counter = 0

    for item01, item02 in itertools.izip(fq_open(fastq_file01), fq_open(fastq_file02)):
        restriction_enzyme_index = int(umi_len) + int(barcode_len)
        counter += 1

        # Restrict to reads with UMI tag and Sample Barcode of full length and with the restriction site as expected
        if item01[1][restriction_enzyme_index:(restriction_enzyme_index + rs_len)] == restriction_enzyme:
            xtag = fq_get_bcmolid(item01, item02, umi_len, restriction_enzyme_index, rs_len, molid_len)
            header01 = item01[0].rstrip()
            header02 = item02[0].rstrip()

            # Add Sample barcode and Molecular_ID to the header of each read in both fastq files,
            # while we remove the UMI tag and Sample Barcode (as well asd their scores) from the first fastq file
            item01[0] = '%s %s\n' % (header01, xtag) 
            item01[1] = item01[1][restriction_enzyme_index:]
            item01[3] = item01[3][restriction_enzyme_index:]
            item02[0] = '%s %s\n' % (header02, xtag)

            for j in range(4):
                fw.write(item01[j])
                bw.write(item02[j])
        else:
            for j in range(4):
                non_fw.write(item01[j])
                non_bw.write(item02[j])
    fw.close()
    bw.close()
    non_fw.close()
    non_bw.close()

    print('Total reads %s' % counter, file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Add barcode and umi:seq to header')

    parser.add_argument('--fastq_file01', help='fastq file with Reads 01', required=True)
    parser.add_argument('--fastq_file02', help='fastq file with Reads 02', required=True)
    parser.add_argument('--umi_len', help='number of bases in the UMI tags', default=8)
    parser.add_argument('--barcode_len', help='number of bases in the barcode samples', default=8)
    parser.add_argument('--restriction_enzyme', help='restriction enzyme site', required=True)
    parser.add_argument('--molid_len', help='length of genomic sequence following the restriction site that will be used in the tag', default=6)
    parser.add_argument('--out_name01', help='basename for output fastq associated to Read 1, including the ".fastq" file extension', required=True)
    parser.add_argument('--out_name02', help='basename for output fastq associated to Read 2, including the ".fastq" file extension', required=True)
    parser.add_argument('--out_dir', help='directory where to save output files', default='.')

    args = parser.parse_args()

    print("Parameters: number of bases in \n barcode = %s \n UMI tag %s \n molecular id = %s" % (args.barcode_len, args.umi_len, args.molid_len), file=sys.stderr)

    fq_umitag(args.fastq_file01, args.fastq_file02, args.umi_len, args.barcode_len, args.restriction_enzyme, args.molid_len, args.out_name01, args.out_name02, args.out_dir)


if __name__ == "__main__":
    main()
