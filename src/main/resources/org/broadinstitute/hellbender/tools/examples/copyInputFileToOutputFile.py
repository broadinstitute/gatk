# Test script that accepts two input arguments that are file paths.
# Copies the contents of th first file to the second file:
import sys
with open(sys.argv[1]) as fin:
    lines = fin.readlines()
    with open(sys.argv[2], "w") as fout:
        fout.writelines(lines)

