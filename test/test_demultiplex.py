"""
test_demultiplex.py
-------------------------------

Tests for the `demultiplex.py` module

"""

import unittest
import os
import sys
import shutil
import utils

# Include the parent directory in the PYTHONPATH for relative imports
sys.path.append('..')
import demultiplex


TEST_SAMPLE_BARCODES = {'AGGCATGAGATCGC': 'mysample', 'GACTCCTGCGATAT': 'sample2'}
TEST_DATA_FILES = {'read1': 'data/undemultiplexed/undemux.r1.fastq',
                  'read2': 'data/undemultiplexed/undemux.r2.fastq',
                  'index1': 'data/undemultiplexed/undemux.i1.fastq',
                  'index2': 'data/undemultiplexed/undemux.i2.fastq'}
TEST_OUTPUT_PATH = 'output'
TEST_MIN_READS = 1000
CORRECT_DEMULTIPLEX_OUTPUT_FOLDER = 'data/demultiplexed'

class TestDemultiplex(unittest.TestCase):

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

        self.assertTrue(utils.checkFolderEquality(TEST_OUTPUT_PATH, CORRECT_DEMULTIPLEX_OUTPUT_FOLDER))


    def tearDown(self):
        # Delete the output folder and the results
        shutil.rmtree(TEST_OUTPUT_PATH)

if __name__ == '__main__':
    unittest.main()
