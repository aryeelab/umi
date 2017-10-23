"""
test_demultiplex.py
-------------------------------

Tests for the `demult_ct.py` module for demultiplixing circulating tumor
samples by sample name and barcode

"""

import unittest
import os
import sys
import shutil
import utils

# Include the parent directory in the PYTHONPATH for relative imports
sys.path.append('..')
import demultiplex_nobuff as demultiplex


P5_SAMPLE_BARCODES = {'GCGATAT': 'P51', 'AGATCGC': 'P52'}
P7_SAMPLE_BARCODES = {'AGGCATG': 'P71', 'GACTCCT': 'P72'}
TEST_DATA_FILES = {'read1': 'test/data/undemultiplexed/undemux.r1.fastq',
                  'read2': 'test/data/undemultiplexed/undemux.r2.fastq',
                  'index1': 'test/data/undemultiplexed/undemux.i1.fastq',
                  'index2': 'test/data/undemultiplexed/undemux.i2.fastq'}
TEST_OUTPUT_PATH = 'output'
TEST_MIN_READS = 1000
CORRECT_DEMULTIPLEX_OUTPUT_FOLDER = 'test/data/demultiplex_nobuff'

class TestDemultiplex(unittest.TestCase):

    def setUp(self):
        # Create the output folder
        if not os.path.exists(TEST_OUTPUT_PATH):
            os.makedirs(TEST_OUTPUT_PATH)


    def testDemultiplexTestCase(self):
        # Run the demultiplex module on the test data
        demultiplex.demultiplex(TEST_DATA_FILES['read1'],
                                TEST_DATA_FILES['read2'],
                                TEST_DATA_FILES['index1'],
                                TEST_DATA_FILES['index2'],
                                P5_SAMPLE_BARCODES,
                                P7_SAMPLE_BARCODES,
                                TEST_OUTPUT_PATH,
                                min_reads=TEST_MIN_READS)

        self.assertTrue(utils.checkFolderEquality(TEST_OUTPUT_PATH, CORRECT_DEMULTIPLEX_OUTPUT_FOLDER))


    def tearDown(self):
        # Delete the output folder and the results
        shutil.rmtree(TEST_OUTPUT_PATH)

if __name__ == '__main__':
    unittest.main()
