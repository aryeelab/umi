"""
test_consolidate.py
-------------------------------

Tests for the `consolidate.py` module

"""

import unittest
import os
import sys
import shutil
import utils

# Include the parent directory in the PYTHONPATH for relative imports
sys.path.append('..')
import consolidate

TEST_OUTPUT_PATH = 'output'
TEST_DATA_FILES = {'read1': 'data/umitagged/mysample.r1.umitagged.fastq',
                  'read2': 'data/umitagged/mysample.r2.umitagged.fastq',
                  'read1_out': os.path.join(TEST_OUTPUT_PATH, 'mysample.r1.consolidated.fastq'),
                  'read2_out': os.path.join(TEST_OUTPUT_PATH, 'mysample.r2.consolidated.fastq')}
TEST_MIN_QUAL = 15
TEST_MIN_FREQ = 0.9
CORRECT_UMITAGGED_OUTPUT_FOLDER = 'data/consolidated'

class TestConsolidate(unittest.TestCase):

    def setUp(self):
        # Create the output folder
        os.makedirs(TEST_OUTPUT_PATH)


    def testConsolidateTestCase(self):
        # Run the consolidation module on the test data
        consolidate.consolidate(TEST_DATA_FILES['read1'],
                                TEST_DATA_FILES['read1_out'],
                                TEST_MIN_QUAL,
                                TEST_MIN_FREQ)

        consolidate.consolidate(TEST_DATA_FILES['read2'],
                                TEST_DATA_FILES['read2_out'],
                                TEST_MIN_QUAL,
                                TEST_MIN_FREQ)

        self.assertTrue(utils.checkFolderEquality(TEST_OUTPUT_PATH, CORRECT_UMITAGGED_OUTPUT_FOLDER))


    def tearDown(self):
        # Delete the output folder and the results
        shutil.rmtree(TEST_OUTPUT_PATH)

if __name__ == '__main__':
    unittest.main()
