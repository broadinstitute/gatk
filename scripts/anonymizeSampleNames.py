#!/usr/bin/python

# Tool to remove sample names from a vcf file and replace with generic Sample_1 Sample2
# It removes ##GATKCommandLine headerlines and replaces the sample names in the vcf header
# Always check files manually to be sure there aren't any other sources of personally identifiable data

import sys

if len(sys.argv) is not 3:
    exit("Usage:\n\tAnonymizeSampleNames <input.vcf> <output.vcf>\n\nInput was:\n\t" + " ".join(sys.argv[1:]))

input = sys.argv[1]
output = sys.argv[2]


def fixSampleNames(line):
    if line.startswith("#CHROM"):
        fields = line.split("\t")
        samples = fields[9:]
        for (i, name) in enumerate(samples):
            fields[9 + i] = "Sample_" + str(i)
        return "\t".join(fields) + "\n"
    elif line.startswith("##GATKCommandLine"):
        return ""
    else:
        return line


with open(input, 'r') as i, open(output, 'w') as o:
    for line in i:
        o.write(fixSampleNames(line))
