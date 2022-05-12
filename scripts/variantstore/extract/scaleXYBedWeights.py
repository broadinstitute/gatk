"""
Script to scale X and Y chromosome weights in a BED file to allow for more balanced extract times.

Usage:

cat original.bed | python3 scaleXYBedWeights.py 4 > scaledXY.bed

"""
import sys


def scaleXY(multiplier):
    while True:
        line = sys.stdin.readline()
        if not line:
            break

        line = line.rstrip('\n')

        if line.startswith('chrX') or line.startswith('chrY'):
            fields = line.split('\t')
            weight = fields[-1]
            fields[-1] = str(int(weight) * multiplier)
            print('\t'.join(fields))
        else:
            print(line)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Must supply a single argument with the multiplier size for X and Y chromosomes")
        exit(1)

    multiplier=sys.argv[1]
    try:
        multiplier=int(sys.argv[1])
        if multiplier < 1:
            print("Multiplier must be a positive integer")
            exit(1)
    except ValueError:
        print("Multiplier value must be an integer")
        exit(1)

    scaleXY(multiplier)
