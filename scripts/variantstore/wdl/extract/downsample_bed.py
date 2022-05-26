from collections import defaultdict
import sys


counts = defaultdict(lambda: 0)

for line in map(str.rstrip, sys.stdin):
    xsome = line.split('\t')[0]

    if counts[xsome] < 5:
        counts[xsome] = counts[xsome] + 1
        print(line)
