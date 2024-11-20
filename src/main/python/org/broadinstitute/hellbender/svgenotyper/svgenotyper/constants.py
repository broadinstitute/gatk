from enum import Enum

INPUT_GENOTYPE_FIELDS = ['PE', 'SSR', 'ESR', 'NCN']

MAX_PE_COUNT = 50
MAX_SR_COUNT = 50
MAX_GT_LOD = 99

# Ploidy of mean depth table input
DEPTH_PLOIDY = 2


class SVTypes(Enum):
    DEL = 0
    DUP = 1
    INS = 2
    INV = 3
    BND = 4
