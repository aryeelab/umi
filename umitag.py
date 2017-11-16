from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse
import subprocess
import multiprocessing as mp
import glob

__author__ = 'Martin Aryee, Allison MacLeay'

# python ~/work/projects/umi/umitag.py --read1_in bcl/Undetermined_S0_L001_R1_001.fastq.gz --read2_in bcl/Undetermined_S0_L001_R2_001.fastq.gz --read1_out umi1.fastq.gz --read2_out umi2.fastq.gz --index1 bcl/Undetermined_S0_L001_I1_001.fastq.gz --index2 bcl/Undetermined_S0_L001_I2_001.fastq.gz

# args = {'out_dir':'/PHShome/ma695/tmp', 'min_reads':10}
# base = '/data/joung/sequencing_bcl/131007_M01326_0075_000000000-A6B33/Data/Intensities/BaseCalls'
# args['read1_in'] = os.path.join(base, 'Undetermined_S0_L001_R1_001.fastq.gz')
# args['read2_in'] = os.path.join(base, 'Undetermined_S0_L001_R2_001.fastq.gz')
# args['read1_out'] = os.path.join(out_dir, 'Undetermined_S0_L001_R1_001.umi.fastq.gz')
# args['read2_out'] = os.path.join(out_dir, 'Undetermined_S0_L001_R2_001.umi.fastq.gz')
# args['index1'] = os.path.join(base, 'Undetermined_S0_L001_I1_001.fastq.gz')
# args['index2'] = os.path.join(base, 'Undetermined_S0_L001_I2_001.fastq.gz')


def fq(file, start, stop):
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')
    ct = start / 4
    with fastq as f:
        for _ in range(start):
            f.readline()
        while (ct * 4) < stop - 1:
            ct += 1
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


def process_fq(read1_out, read2_out, read1, read2, index1, index2, pattern, start, stop):
    r1_umitagged_unsorted_file = '{}_{}.tmp'.format(read1_out, start)
    r2_umitagged_unsorted_file = '{}_{}.tmp'.format(read2_out, start)

    # Create UMI-tagged R1 and R2 FASTQs
    r1_umitagged = open(r1_umitagged_unsorted_file, 'w')
    r2_umitagged = open(r2_umitagged_unsorted_file, 'w')
    for r1, r2, i1, i2 in itertools.izip(fq(read1, start, stop), fq(read2, start, stop), fq(index1, start, stop), fq(index2, start, stop)):
        # Create molecular ID by concatenating molecular barcode and beginning of r1 read sequence
        molecular_id, r1[1], r2[1], r1[3], r2[3] = get_molecular_barcode(r1, r2, i1, i2, pattern)
        # Add molecular id to read headers
        r1[0] = '%s %s\n' % (r1[0].rstrip(), molecular_id)
        r2[0] = '%s %s\n' % (r2[0].rstrip(), molecular_id)
        r1[0] = r1[0].replace(' ', ':')
        r2[0] = r2[0].replace(' ', ':')
        for line in r1:
            r1_umitagged.write(line)
        for line in r2:
            r2_umitagged.write(line)
    r1_umitagged.close()
    r2_umitagged.close()
    return r1_umitagged_unsorted_file, r2_umitagged_unsorted_file

def get_numlines(fpath):
    ct = 0
    with open(fpath, 'rb') as fh:
        for _ in fh:
            ct += 1
    return ct

def merge_output(res, num_procs):
    r1_umitagged_unsorted_file = None
    r2_umitagged_unsorted_file = None
    for r in res:
        r1, r2 = r.get()
        if r1_umitagged_unsorted_file is None:
            r1_umitagged_unsorted_file = r1
            r2_umitagged_unsorted_file = r2
        else:
            cmds = ['cat {} >> {}'.format(r1, r1_umitagged_unsorted_file),
                    'cat {} >> {}'.format(r2, r2_umitagged_unsorted_file)]
            if num_procs > 1:
                pool = mp.Pool(processes=2)
                procs = [pool.apply_async(os.system, args=(arg, )) for arg in cmds]
                for proc in procs:
                    proc.get()
            else:
                os.system(cmds[0])
                os.system(cmds[1])
    return r1_umitagged_unsorted_file, r2_umitagged_unsorted_file


def umitag(read1, read2, index1, index2, read1_out, read2_out, out_dir, pattern, num_procs):

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    read1_out = os.path.join(out_dir, os.path.basename(read1_out))
    read2_out = os.path.join(out_dir, os.path.basename(read2_out))

    if not index1:  # placeholder
        index1 = read1
    if not index2:  # placeholder
        index2 = read2
    num_lines = get_numlines(index1)
    chunk_size = num_lines / num_procs
    if num_lines % num_procs != 0:  # math.ceil()
        chunk_size += 1
    diff = chunk_size % 4
    chunk_size += (4 - diff)
    pool = mp.Pool(processes=num_procs)
    res = [pool.apply_async(process_fq, args=(read1_out, read2_out, read1, read2, index1, index2, pattern, chunk * chunk_size, (chunk + 1) * chunk_size - 1)) for chunk in range(num_procs)]
    r1_umitagged_unsorted_file, r2_umitagged_unsorted_file = merge_output(res, num_procs)
    # Sort fastqs based on molecular barcode
    cmd1 = 'cat ' + r1_umitagged_unsorted_file + ' | paste - - - - | sort -k1,1 | tr "\t" "\n" >' + read1_out
    cmd2 = 'cat ' + r2_umitagged_unsorted_file + ' | paste - - - - | sort -k1,1 | tr "\t" "\n" >' + read2_out
    if num_procs > 1:
        pool = mp.Pool(processes=2)
        procs = [pool.apply_async(subprocess.check_call, args=(cmd, ), kwds=dict(shell=True, env=os.environ.copy())) for cmd in [cmd1, cmd2]]
    else:
        subprocess.check_call(cmd1, shell=True, env=os.environ.copy())
        subprocess.check_call(cmd2, shell=True, env=os.environ.copy())
    for proc in procs:
        proc.get()
    tmp_files = glob.glob(os.path.join(os.path.dirname(r1_umitagged_unsorted_file), '*.tmp'))
    for tmp_file in tmp_files:
        os.remove(tmp_file)


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
    parser.add_argument('--threads', default=1)
    args = vars(parser.parse_args())

    r1_pat, r2_pat, i1_pat, i2_pat = args['pattern'].split(',')
    if (i1_pat != '' or i2_pat != '') and (None in [args['index1'], args['index2']]):
        print('Index files are required when pattern is defined to use indexes!')

    umitag(args['read1_in'], args['read2_in'], args['index1'], args['index2'],
           args['read1_out'], args['read2_out'], args['out_dir'], args['pattern'], int(args['threads']))

if __name__ == '__main__':
    main()
