"""
test_consolidate.py
-------------------------------

Tests for the `consolidate.py` module

"""

import unittest
import os
import sys
import inspect
import shutil
import filecmp

# Include the parent directory in the PYTHONPATH for relative imports
sys.path.append('..')
import consolidate


TEST_SAMPLE_BARCODES = {'AGGCATGAGATCGC': 'mysample', 'GACTCCTGCGATAT': 'sample2'}
TEST_DATA_FILES = {'read1': 'data/undemux/undemux.r1.fastq',
                  'read2': 'data/undemux/undemux.r2.fastq',
                  'index1': 'data/undemux/undemux.i1.fastq',
                  'index2': 'data/undemux/undemux.i2.fastq'}
TEST_OUTPUT_PATH = 'output'
TEST_MIN_READS = 1000
CORRECT_DEMULTIPLEX_OUTPUT_FOLDER = 'data/undemux/demux_results'

class TestConsolidate(unittest.TestCase):

    def setUp(self):
        # Create the output folder
        os.makedirs(TEST_OUTPUT_PATH)


    def testDemultiplexTestCase(self):
        # Run the demultiplex module on the test data
        demultiplex.demultiplex(TEST_DATA_FILES['read1'],
                                TEST_DATA_FILES['read2'],
                                TEST_DATA_FILES['index1'],
                                TEST_DATA_FILES['index2'],
                                TEST_SAMPLE_BARCODES,
                                TEST_OUTPUT_PATH,
                                min_reads=TEST_MIN_READS)

        self.assertTrue(checkFolderEquality(TEST_OUTPUT_PATH, CORRECT_DEMULTIPLEX_OUTPUT_FOLDER))


    def tearDown(self):
        # Delete the output folder and the results
        shutil.rmtree(TEST_OUTPUT_PATH)


def checkFolderEquality(folder1, folder2):
    """
    Given two folders, check if there are the same number of files,
    that the names of files are the same, and that the files with the same
    names are the same.
    """

    folder1_files = os.listdir(folder1)
    folder2_files = os.listdir(folder2)

    if set(folder1_files) != set(folder2_files):
        return False

    for f in folder1_files:
        if not filecmp.cmp(os.path.join(folder1, f), os.path.join(folder2, f)):
            return False

    return True

if __name__ == '__main__':
    unittest.main()
