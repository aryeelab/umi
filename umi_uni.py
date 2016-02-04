"""
umi_uni.py
"""

from __future__ import print_function

___author___ = 'Martin Aryee'

### modules
import gzip
import re
import sys
import os
import argparse

""" retrieve the file, allowing .gz formats """
def fq_open(arch):
    if re.search(".gz$", arch):
        return( gzip.open(arch, "rb") )
    else:
        return( open(arch, "r") )

""" yielding the 4 lines of each reading as a generator object """
def fq_yield(arch):
    openarch = fq_open(arch)
    while True:
        line1 = openarch.readline()
        if not line1:
            break
        line2 = openarch.readline()
        line3 = openarch.readline()
        line4 = openarch.readline()
        yield[line1, line2, line3, line4]

""" retrieve barcode and molecular_id of a "yielded" fastq element as two strings separated by an empty space"""
def fq_getbcumi(element, cw, cw_index, barcode_len, id_len):
    cutindex = element[1].find(cw)
    barcode = element[1][cutindex-int(barcode_len):cutindex]
    id_pos = cutindex+len(cw)
    #umi_(first id_len base pairs after restriction enzyme in the actual genomic sequence
    molecular_id = element[1][:cutindex-int(barcode_len)] +':'+ element[1][id_pos:id_pos+int(id_len)] 
    out = '%s %s' % (barcode, molecular_id)
    return(out)

""" re-tagged first line in fastq file: seq_id line_id barcode umi
      score of barcode and umi are removed from line 4"""
def fq_umitag(fastq_file01, cw, cw_index, barcode_len, id_len, out_name01, out_dir):

    if int(cw_index) - int(barcode_len) <= 0:
        print('ERROR: dimensions of barcode and restriction enzyme are not compatible.', file=sys.stderr)
        return()
    if int(barcode_len)<=0 or int(id_len)<=0 or int(cw_index)<=0:
        print('ERROR: numerical input values can not be negative.', file=sys.stderr)
        return()        
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out_name01_path = os.path.join(out_dir, out_name01)
    mout_name01_path = os.path.join(out_dir, 'NONstored_' + out_name01)
    
    fw = open(out_name01_path, "w")
    mfw = open(mout_name01_path, "w")
    bc_list = []

    for item01 in fq_yield(fastq_file01):
        cwindex = item01[1].find(cw)
        
        # condition necessary to avoid missing cutword and length of  barcode-umi = index position of cw
        if cwindex == int(cw_index):
            xtag = fq_getbcumi(item01, cw, cw_index, barcode_len, id_len) 
            header01 = item01[0].strip()
            barcode =  xtag.split()[0]
            
            if barcode not in bc_list:
                bc_list = bc_list + [barcode]
            
            item01[0] = '%s %s\n' % (header01, xtag) 
            item01[1] = item01[1][cwindex:]
            item01[3] = item01[3][cwindex:]
            for j in range(4):
                fw.write(item01[j])
        else:
            for j in range(4):
                mfw.write(item01[j])

    fw.close()
    mfw.close()


def main():
    parser = argparse.ArgumentParser(description='Add barcode_id and umi_tag to header.')
    parser.add_argument('--fastq_file', help='original fastq file', required=True)
    parser.add_argument('--cw', help='restriction enzyme', required=True)
    parser.add_argument('--cw_index', help='zero-base position of restriction enzyme in the base par sequence', default=16)
    parser.add_argument('--barcode_len', help='number of bases in the barcodes', default=8)
    parser.add_argument('--id_len', help='length of genenomic sequence following the restriction enzyme that to will be used as part of the umi tag', default=6)
    parser.add_argument('--out_name', help='output fastq name', required=True)
    parser.add_argument('--out_dir', help='directory where to save output files', default='.')
    args = parser.parse_args()

    print("Parameters: \n number of bases in each barcode = %s \n restriction enzyme (0-based) starting position = %s \n number of bases to be taken for the molecular id = %s" \
          % (args.barcode_len, args.cw_index, args.id_len), file=sys.stderr)
    fq_umitag(args.fastq_file, args.cw, args.cw_index, args.barcode_len , args.id_len, args.out_name, args.out_dir)

if __name__ == "__main__":
    main()
