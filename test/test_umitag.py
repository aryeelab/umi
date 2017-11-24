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
import umitag

TEST_OUTPUT_PATH = 'output'
TEST_DATA_FILES = {'read1': 'data/demultiplexed/mysample.r1.fastq',
                  'read2': 'data/demultiplexed/mysample.r2.fastq',
                  'index1': 'data/demultiplexed/mysample.i1.fastq',
                  'index2': 'data/demultiplexed/mysample.i2.fastq',
                  'read1_out': os.path.join(TEST_OUTPUT_PATH, 'mysample.r1.umitagged.fastq'),
                  'read2_out': os.path.join(TEST_OUTPUT_PATH, 'mysample.r2.umitagged.fastq')}
CORRECT_UMITAGGED_OUTPUT_FOLDER = 'data/umitagged'

class TestUMITag(unittest.TestCase):

    def setUp(self):
        """ Create the output folder """
        if not os.path.exists(TEST_OUTPUT_PATH):
           os.makedirs(TEST_OUTPUT_PATH)
        with open(os.path.join(TEST_OUTPUT_PATH, 'diffsample.r1.umitagged.fastq'), 'w') as fh:
           fh.write(' ')

    def testSortVersion(self):
        """ try to run sort --version """
        assert(umitag.get_sort_opts() in ['', ' -V '])

    def testUMITagTestCase(self):
        """ Run the umitag module on the test data """
        print(os.getcwd())
        umitag.umitag(TEST_DATA_FILES['read1'],
                                TEST_DATA_FILES['read2'],
                                TEST_DATA_FILES['index1'],
                                TEST_DATA_FILES['index2'],
                                TEST_DATA_FILES['read1_out'],
                                TEST_DATA_FILES['read2_out'],
                                TEST_OUTPUT_PATH, '8B12H,,,', 8)

        self.assertTrue(utils.checkFolderEquality(TEST_OUTPUT_PATH, CORRECT_UMITAGGED_OUTPUT_FOLDER))

    def tearDown(self):
        """ Delete the output folder and the results """
        shutil.rmtree(TEST_OUTPUT_PATH)


if __name__ == '__main__':
    unittest.main()
