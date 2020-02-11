from enum import Enum

GENOTYPE_FIELDS = ['PE', 'SSR', 'ESR', 'NCN']

MAX_PE_COUNT = 50
MAX_SR_COUNT = 50


class SVTypes(Enum):
    DEL = 0
    DUP = 1
    INS = 2
    INV = 3
    BND = 4
