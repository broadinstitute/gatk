from argparse import ArgumentParser, RawDescriptionHelpFormatter
from copy import deepcopy
import os
import subprocess
import pysam
import random

# DO NOT SEED
#random.seed("3337")

def parse_options(description, epilog):
    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument('input_file', type=str, help='')
    parser.add_argument('contigs', type=str, help='The contigs to PRESERVE separated by "," ... THESE MUST BE SORTED (e.g. 5,20,21 is okay ; 20,5,21 is not okay)')
    parser.add_argument('output_file', type=str, help='')

    # Process arguments
    args = parser.parse_args()

    return args


def main():
    args = parse_options("Tool to reheader a (copy of a) bam to only the specified contigs.  EXPERIMENTAL and UNSUPPORTED", """WARNING:  No error checking.
    Assumes that the bam is coordinate sorted and paired-end.
    java -jar picard.jar ReplaceSamHeader HEADER=tmp_header.sam I=<(cat <(head -n1  tmp_header.sam ) <(samtools view tumor_1_foo.bam)) O=yossi.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

    THIS SCRIPT CAN BE VERY SLOW ON LARGE BAMS

    This script is mostly meant for generating test bams that validate (even in Picard).

    Example: python reheader_bam.py some.bam 20,21 some_20_21_only.bam""")

    input_file = args.input_file
    output_file = args.output_file
    contigs = args.contigs.split(",")
    contigs = set([c.strip() for c in contigs])

    # 1) BAM -> Filtered SAM w/ old header
    random_num = random.random()
    samfile_out_filename_step1 = os.path.dirname(os.path.abspath(output_file)) + "/tmp_step_1_" + str(random_num) + ".sam"

    bamfile_in = pysam.AlignmentFile(input_file, 'rb')

    samfile_out = pysam.AlignmentFile(samfile_out_filename_step1, 'wh', header=bamfile_in.header)
    for i,r in enumerate(bamfile_in):
        if r.next_reference_name in contigs:
            try:
                s_mate = bamfile_in.mate(r)
                samfile_out.write(r)
            except ValueError as ve:
                if ve.message.find("fetch called") != -1 or ve.message.find("fetching by") != -1:
                    raise ve

        if i % 100 == 0:
            print(str(i) + " reads... (" + input_file + ")")

    samfile_out.close()

    # 2) SAM Filtered -> new header only SAM file
    # samfile_out_filename is the output from previous step which is input now.
    samfile_in = pysam.AlignmentFile(samfile_out_filename_step1, 'r')
    samfile_out_filename_step2 = os.path.dirname(os.path.abspath(output_file)) + "/tmp_step_2_" + str(random_num) + ".sam"

    # Cosntruct the header-only sam file
    old_header = deepcopy(samfile_in.header)
    new_sq = [c for c in old_header['SQ'] if c['SN'] in contigs]
    new_header = deepcopy(old_header)
    new_header['SQ'] = new_sq

    samfile_out = pysam.AlignmentFile(samfile_out_filename_step2, 'wh', header=new_header)
    samfile_out.close()


    # 3) SAM Filtered + new header only SAM file -> BAM
    bamfile_out_filename_step3 = output_file
    header_only_in = pysam.AlignmentFile(samfile_out_filename_step2, 'r')

    # Create a dictionary, so that we can lookup the SN and get the index.
    sn_to_index_dict = dict()
    for i, entry in enumerate(header_only_in.header['SQ']):
        sn_to_index_dict[entry['SN']] = i

    samfile_in = pysam.AlignmentFile(samfile_out_filename_step1, 'r')
    bamfile_out = pysam.AlignmentFile(bamfile_out_filename_step3, 'wb', header=header_only_in.header)
    for r in samfile_in:
        r.reference_id = sn_to_index_dict[r.reference_name]
        r.next_reference_id = sn_to_index_dict[r.next_reference_name]
        bamfile_out.write(r)

    bamfile_out.close()

    subprocess.call("samtools index " + bamfile_out_filename_step3, shell=True)

if __name__ == "__main__":
    main()
