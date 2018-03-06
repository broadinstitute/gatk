# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Definitions ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import re

TENSOR_MAPS_2D = ['read_tensor']
TENSOR_MAPS_1D = ['reference']
TENSOR_SUFFIX = '.hd5'

DNA_SYMBOLS = {'A':0, 'C':1, 'G':2, 'T':3}
INPUTS_INDEL = {'A':0, 'C':1, 'G':2, 'T':3, '*':4}
INPUT_SYMBOLS = {
    'dna' : DNA_SYMBOLS,
    'dna_indel' : INPUTS_INDEL,
}

# Base calling ambiguities, See https://www.bioinformatics.org/sms/iupac.html
AMBIGUITY_CODES = {
    'K':[0, 0, 0.5, 0.5], 'M':[0.5, 0.5, 0, 0], 'R':[0.5, 0, 0, 0.5], 'Y':[0, 0.5, 0.5, 0], 'S':[0, 0.5, 0, 0.5],
    'W':[0.5, 0, 0.5, 0], 'B':[0,0.333,0.333,0.334], 'V':[0.333,0.333,0,0.334], 'H':[0.333,0.333,0.334,0],
    'D':[0.333,0,0.333,0.334], 'X':[0.25,0.25,0.25,0.25], 'N':[0.25,0.25,0.25,0.25]
}


# Annotation sets
ANNOTATIONS = {
    '_' : [], # Allow command line to unset annotations
    'best_practices_w_qual' : ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'QUAL', 'ReadPosRankSum'],
    'best_practices' : ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'ReadPosRankSum'],
    'annotations' : ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'ReadPosRankSum'],
    'm2':['AF', 'AD_0', 'AD_1', 'MBQ', 'MFRL_0', 'MFRL_1', 'MMQ', 'MPOS'],
    'combine': ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'ReadPosRankSum', 'AF', 'AD_0', 'AD_1', 'MBQ', 'MFRL_0', 'MFRL_1', 'MMQ', 'MPOS'],
    'gnomad': ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'ReadPosRankSum', 'DP_MEDIAN', 'DREF_MEDIAN', 'GQ_MEDIAN', 'AB_MEDIAN'],
}

SNP_INDEL_LABELS = {'NOT_SNP':0, 'NOT_INDEL':1, 'SNP':2, 'INDEL':3}

CODE2CIGAR = 'MIDNSHP=XB'
CIGAR2CODE = dict([y, x] for x, y in enumerate(CODE2CIGAR))
CIGAR_CODE = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4}
CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")

SKIP_CHAR = '~'
INDEL_CHAR = '*'
SEPARATOR_CHAR = '\t'

MAPPING_QUALITY_MAX = 60.0 # Mapping qualities from BWA are typically capped at 60
READ_FLAGS = 12 # Total number of read flags, actual flags used is determined by the tensor map

