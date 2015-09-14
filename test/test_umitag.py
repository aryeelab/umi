"""
test_umitag.py
-------------------------------

Tests for the `umitag.py` module

"""

import unittest
import os
import sys
import shutil
import utils

# Include the parent directory in the PYTHONPATH for relative imports
sys.path.append('..')
import demultiplex


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
