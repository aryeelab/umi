import os
import sys
import inspect
import filecmp

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