from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse
import subprocess
import sys

__author__ = 'Martin Aryee, Allison MacLeay'


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


def get_molecular_barcode(read1, read2, index1, index2, strategy='8B12X,,,'):
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
            elif btype == 'X':  # Handle - skip this
                pass
            pos += num
        if i == 0:
            r1_out = read[1][pos:]
            q1_out = read[3][pos:]
        elif i == 1:
            r2_out = read[1][pos:]
            q2_out = read[3][pos:]
    return barcode, r1_out, r2_out, q1_out, q2_out


def tag_query_name(query_name, molecular_id):
    space = query_name.split(' ')
    pts = space[0].split(':')
    pts[2] = '{}_{}'.format(pts[2], molecular_id)
    query_name = '{}\n'.format(':'.join(pts))
    if len(space) > 1:
        query_name = '{} {}'.format(query_name.replace('\n', ''), ' '.join(space[1:]))
    return query_name


def process_fq(read1_out, read2_out, read1, read2, index1, index2, pattern):
    r1_umitagged_unsorted_file = '{}.tmp'.format(read1_out)
    r2_umitagged_unsorted_file = '{}.tmp'.format(read2_out)

    # Create UMI-tagged R1 and R2 FASTQs
    r1_umitagged = open(r1_umitagged_unsorted_file, 'w')
    r2_umitagged = open(r2_umitagged_unsorted_file, 'w')
    try:
        # support python 3
        zip_func = itertools.izip if sys.version_info[0] < 3 else zip
        for r1, r2, i1, i2 in zip_func(fq(read1), fq(read2), fq(index1), fq(index2)):
            # Create molecular ID by concatenating molecular barcode and beginning of r1 read sequence
            molecular_id, r1[1], r2[1], r1[3], r2[3] = get_molecular_barcode(r1, r2, i1, i2, pattern)
            # Add molecular id to read headers
            r1[0] = tag_query_name(r1[0], molecular_id)
            r2[0] = tag_query_name(r2[0], molecular_id)
            for line in r1:
                r1_umitagged.write(line)
            for line in r2:
                r2_umitagged.write(line)
    except Exception as e:
        print(e)
        raise e
    return r1_umitagged_unsorted_file, r2_umitagged_unsorted_file


def sort_fastqs(r_umitagged_unsorted_file, sort_opts, read_out):
    """ sort by query name.
        https://edwards.sdsu.edu/research/sorting-fastq-files-by-their-sequence-identifiers/
    """
    cmd = 'cat {0} | paste - - - - | sort -k1,1 {1} | tr "\t" "\n" > {2}; rm {0}'.format(r_umitagged_unsorted_file, sort_opts,
                                                                              read_out)
    sort_output = subprocess.check_output(cmd, shell=True, env=os.environ.copy())
    return sort_output


def get_sort_opts():
    try:
        version = float(re.match('sort \(\w+ \w+\) ([\d\.]*)',
                                 subprocess.check_output('sort --version', shell=True)).group(1))
        if version > 5.93:
            return ' -V '
    except Exception as e:
        print(e)
    return ''


def umitag(read1, read2, index1, index2, read1_out, read2_out, out_dir, pattern):

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    read1_out = os.path.join(out_dir, os.path.basename(read1_out))
    read2_out = os.path.join(out_dir, os.path.basename(read2_out))

    r1_umitagged_unsorted_file = read1_out + '.tmp'
    r2_umitagged_unsorted_file = read2_out + '.tmp'

    if not index1:  # placeholder
        index1 = read1
    if not index2:  # placeholder
        index2 = read2

    r1_umitagged_unsorted_file, r2_umitagged_unsorted_file = process_fq(r1_umitagged_unsorted_file,
                                                                        r2_umitagged_unsorted_file,
                                                                        read1, read2, index1, index2, pattern)
    sort_opts = get_sort_opts()
    args_list = [(r1_umitagged_unsorted_file, sort_opts, read1_out), (r2_umitagged_unsorted_file, sort_opts, read2_out)]

    # Sort fastqs based on molecular barcode
    for args in args_list:
        sort_fastqs(*args)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read1_in', required=True)
    parser.add_argument('--read2_in', required=True)
    parser.add_argument('--read1_out', required=True)
    parser.add_argument('--read2_out', required=True)
    parser.add_argument('--index1')
    parser.add_argument('--index2')
    parser.add_argument('--pattern', default='8B12X,,,')
    parser.add_argument('--out_dir', default='.')
    args = vars(parser.parse_args())

    r1_pat, r2_pat, i1_pat, i2_pat = args['pattern'].split(',')
    if (i1_pat != '' or i2_pat != '') and (None in [args['index1'], args['index2']]):
        print('Index files are required when pattern is defined to use indexes!')

    umitag(args['read1_in'], args['read2_in'], args['index1'], args['index2'],
           args['read1_out'], args['read2_out'], args['out_dir'], args['pattern'])

if __name__ == '__main__':
    main()
