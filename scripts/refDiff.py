#!/usr/bin/env python
"""
SYNOPSIS

    %prog [options] REF_FASTA_1 REF_FASTA_2 [REF_FASTA_3 ...]

DESCRIPTION

    Get differences in reference sequences from the same species.
    For instance this will allow you to see what contigs are different between HG19 and b37 (if any)
    with base-level detail where applicable.

    Each reference file must have a sequence dictionary (.dict) file with the same base name in the same
    directory as the FASTA file itself.

EXAMPLES

    %prog ucsc.hg19.fasta Homo_sapiens_assembly19.fasta
    %prog ucsc.hg19.fasta Homo_sapiens_assembly19.fasta GRCh37.p13.genome.fasta

EXIT STATUS

    0 - Success
    1 - Exception raised during execution
    2 - Arguments did not validate

AUTHOR

    Jonn Smith <jonn@broadinstitute.org>

LICENSE

    BSD 3-Clause License.  See LICENSE file for details (top level of repo).
"""

#########################################################################################################
# Style Note: Style for python scripts should follow the Python Software Fonudation's guide found here: #
#             https://www.python.org/dev/peps/pep-0008/                                                 #
#########################################################################################################

################################################################################
# Import section:
####################

# Default imports (for template)
import signal
import sys
import argparse
import re
import time
import logging
import traceback
import os

# Custom imports for script:
import pysam
import csv
from collections import namedtuple
import re

################################################################################
# Built-in Module Vars:
####################

_START_TIME = time.time()
_VERSION = 0.01
_LOG_LEVEL = logging.NOTSET

_PROG_NAME = re.sub('^.*/', '', sys.argv[0])

# Get description info from doc at top:
_epilog_regex = re.compile(r'.*(EXAMPLES.*)', re.MULTILINE | re.DOTALL)
_EPILOG = _epilog_regex.match(globals()['__doc__']).groups(0)[0].replace('%prog', _PROG_NAME)

_description_regex = re.compile(r'.*DESCRIPTION(\s+.*?)EXAMPLES.*', re.MULTILINE | re.DOTALL)
_DESCRIPTION = _description_regex.match(globals()['__doc__']).groups(0)[0]

LOGGER = logging.getLogger(__name__)

########################################
# Signal Handler:
####################


def _interrupt_handler(_, __):
    """Handle the SIGINT signal."""
    LOGGER.warning('Signal Interrupt (CTRL+C) received.  Quitting hard.')
    sys.exit()


signal.signal(signal.SIGINT, _interrupt_handler)


########################################
# Arguments Handling:
####################


def _setup_argument_parser():
    """
    Create the options parser to get at command-line options.
    Returns an OptionParser object ready to parse the given command-line options.
    """
    _parser = argparse.ArgumentParser(description=_DESCRIPTION, epilog=_EPILOG,
                                      formatter_class=argparse.RawTextHelpFormatter)

    log_levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    _parser.add_argument("-v", "--verbosity", action='store', nargs=1, dest="log_level",
                         choices=log_levels, default=['INFO'], required=False,
                         help="Set the logging level.")

    _parser.add_argument('ref_files', metavar='REFERENCE_FASTA', type=str, nargs='+',
                         help="Reference FASTA file to compare.")

    return _parser


def _validate_options_and_args(_args):
    """Perform custom syntactic/semantic validation of options and arguments."""

    valid = True

    if len(_args.ref_files) < 2:
        LOGGER.error("Must pass at least 2 reference FASTA files.  Given: %d", len(_args.ref_files))
        valid = False

    # Make sure our fasta and dict (if required) files exist:
    for fasta_path in _args.ref_files:

        # Check the file first:
        file_valid = True
        if not os.path.exists(fasta_path):
            LOGGER.error("FASTA file does not exist: %s", fasta_path)
            file_valid = False
        elif os.path.isdir(fasta_path):
            LOGGER.error("Given FASTA file path is actually a directory: %s", fasta_path)
            file_valid = False

        # Check the sequence dictionary:
        base_path, _ = os.path.splitext(fasta_path)
        dict_path = base_path + ".dict"
        if not os.path.exists(dict_path):
            LOGGER.error("FASTA sequence dictionary file does not exist: %s", dict_path)
            file_valid = False
        elif os.path.isdir(dict_path):
            LOGGER.error("FASTA sequence dictionary is actually a directory: %s", dict_path)
            file_valid = False

        valid = valid and file_valid

    return valid

########################################
# Custom types:
####################


SeqDictEntry = namedtuple("SeqDictEntry", ['name', 'length', 'md5sum', 'url', 'species'])


########################################
# Helper functions:
####################


def atoi(text):
    """Taken from: https://stackoverflow.com/a/5967539"""
    return int(text) if text.isdigit() else text


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    ---
    Taken from: https://stackoverflow.com/a/5967539
    """
    return [atoi(c) for c in re.split(r'(\d+)', text)]


def get_sequence_dict(dict_path):
    """Return a map containing the names of sequences in the given fasta_path and the attributes of that sequence."""

    seq_dict = dict()

    LOGGER.debug("Reading in FASTA sequence dictionary: %s", dict_path)
    with open(dict_path, 'r') as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            # Ignore headers:
            if row[0] == "@HD":
                continue
            elif row[0] == "@SQ":
                # We have sequence info.  We should store it.

                name = ""
                length = 0
                md5 = ""
                url = ""
                species = ""

                # Not clear that the fields occur in the same order each time:
                for field in row[1:]:
                    if field.startswith("SN"):
                        name = field[3:]
                    elif field.startswith("LN"):
                        length = int(field[3:])
                    elif field.startswith("M5"):
                        md5 = field[3:]
                    elif field.startswith("UR"):
                        url = field[3:]
                    elif field.startswith("SP"):
                        species = field[3:]

                seq_dict[name] = SeqDictEntry(name, length, md5, url, species)

    return seq_dict


def analyze_sequence_dictionaries(_args):
    """Compares the sequence dictionaries of the user-specified FASTA files and returns three dictionaries.
    The first is created such that each item is the base name of the FASTA file and the unique contigs in that FASTA
        FASTA_BASE_NAME -> list( UNIQUE_CONTIG_NAME_1, UNIQUE_CONTIG_NAME_2, UNIQUE_CONTIG_NAME_3, ... )

    The second dictionary will contain a mapping for all contigs that are the same across all given FASTA files:

        FASTA_BASE_NAME -> dict( CONTIG_NAME -> dict( FASTA_BASE_NAME -> IDENTICAL_CONTIG_NAME ) )

    The third dictionary will be a dict of the sequence dictionaries for each given FASTA file.
    """

    seq_dict_map = dict()
    seq_md5sum_set_map = dict()
    seq_md5sum_name_map_map = dict()

    LOGGER.info("Ingesting sequence dictionaries...")

    # Get the sequence dictionaries:
    for fasta_path in _args.ref_files:
        base_path, _ = os.path.splitext(fasta_path)
        base_name = os.path.basename(base_path)

        seq_dict = get_sequence_dict(base_path + ".dict")

        seq_dict_map[base_name] = seq_dict

        # Now get a set of the MD5sums of the sequences so we can find them quickly:
        md5sum_set = {s.md5sum for s in seq_dict.values()}
        seq_md5sum_set_map[base_name] = md5sum_set

        seq_md5sum_name_map_map[base_name] = {s.md5sum: s.name for s in seq_dict.values()}

    LOGGER.debug("Read in dictionary info for %d FASTA files.", len(seq_md5sum_set_map))

    # Now compare them:
    unique_map = dict()
    identical_map = dict()

    LOGGER.info("Analyzing sequence dictionaries...")
    for fasta_base in seq_md5sum_set_map.keys():

        LOGGER.debug("Analyzing %s", fasta_base)
        unique_seq_md5s = set()

        for md5 in seq_md5sum_set_map[fasta_base]:

            identical_contig_map = dict()

            LOGGER.debug("    Contig: %s - %s", seq_md5sum_name_map_map[fasta_base][md5], md5)

            is_unique = True
            for other_fasta_base in seq_md5sum_set_map.keys():
                if fasta_base == other_fasta_base:
                    continue

                if md5 in seq_md5sum_set_map[other_fasta_base]:
                    is_unique = False
                    identical_contig_map[other_fasta_base] = seq_md5sum_name_map_map[other_fasta_base][md5]
                    LOGGER.debug("        Found identical contigs: %s\t%s", seq_md5sum_name_map_map[fasta_base][md5],
                                 seq_md5sum_name_map_map[other_fasta_base][md5])

            if is_unique:
                unique_seq_md5s.add(md5)

            if len(identical_contig_map) != 0:
                if fasta_base not in identical_map:
                    identical_map[fasta_base] = dict()

                contig_name = seq_md5sum_name_map_map[fasta_base][md5]
                identical_map[fasta_base][contig_name] = identical_contig_map

        unique_seqs = [s.name for s in seq_dict_map[fasta_base].values() if s.md5sum in unique_seq_md5s]
        unique_map[fasta_base] = unique_seqs

    return unique_map, identical_map, seq_dict_map


def print_contig_table(unique_contig_map, identical_contig_map, seq_dict_map):
    """Prints a table of the contigs in the given maps according to which are identical or unique
    in each user-supplied reference fasta."""

    col_spacer = '\t'

    # Output should be like this:
    # MD5sum    REF_1    REF_2    [REF_3    REF_4...]

    # Create our format for the output:
    # Maybe someday replace this with an F-string:
    row_format_string = "{:<" + str(32) + "}" + col_spacer

    sorted_ref_names = sorted(seq_dict_map.keys())

    for ref in sorted_ref_names:
        longest_field = len(ref)
        longest_contig_len = max([len(contig_name) for contig_name in seq_dict_map[ref]])
        if longest_field < longest_contig_len:
            longest_field = longest_contig_len

        row_format_string = row_format_string + "{:<" + str(longest_field) + "}" + col_spacer

    # Remove the last spacer:
    row_format_string = row_format_string[:-len(col_spacer)]

    # Print the header:
    header_fields = ["md5sum"] + [name for name in sorted_ref_names]
    print(row_format_string.format(*header_fields))

    # Now print out the identical rows first:

    # We sort first by the order in sorted_ref_names, then by alphabetical order of the contigs of the first reference
    # in sorted_ref_names:
    for ref_name in sorted_ref_names:
        identical_contig_list = sorted(identical_contig_map[ref_name].keys(), key=natural_keys)
        for contig_name in identical_contig_list:
            fields_to_print = [seq_dict_map[ref_name][contig_name].md5sum]
            for ref in sorted_ref_names:
                if ref == ref_name:
                    fields_to_print.append(contig_name)
                elif ref in identical_contig_map[ref_name][contig_name]:
                    fields_to_print.append(identical_contig_map[ref_name][contig_name][ref])
                else:
                    fields_to_print.append("----")
            print(row_format_string.format(*fields_to_print))

    # Now we print out the unique contigs from each ref in the same order:
    for ref_name in sorted_ref_names:
        unique_contig_list = sorted(unique_contig_map[ref_name], key=natural_keys)
        for contig_name in unique_contig_list:
            fields_to_print = [seq_dict_map[ref_name][contig_name].md5sum]
            for ref in sorted_ref_names:
                if ref == ref_name:
                    fields_to_print.append(contig_name)
                else:
                    fields_to_print.append("----")
            print(row_format_string.format(*fields_to_print))

########################################
# Main function:
####################


def _main(_args):
    """Default main implementation."""

    # Do a dictionary comparison:
    unique_contig_map, identical_contig_map, seq_dict_map = analyze_sequence_dictionaries(_args)
    print_contig_table(unique_contig_map, identical_contig_map, seq_dict_map)


################################################################################
# Main section:
####################


if __name__ == '__main__':

    # Parse the args:
    parser = _setup_argument_parser()
    args = parser.parse_args()

    _LOG_LEVEL = getattr(logging, args.log_level[0])
    log_format_string = f"%(asctime)s %(name)s %(levelname)-8s %(message)s"
    logging.basicConfig(level=_LOG_LEVEL, format=log_format_string)

    try:
        LOGGER.info('Invoked with arguments: %s', str(args))

        # Validate our options:
        if not _validate_options_and_args(args):
            LOGGER.critical("Error: options and/or arguments are invalid.")
            sys.exit(2)

        # Do the work here:
        _main(args)

        # Log the time here:
        LOGGER.info('Elapsed time: %2.2f', time.time() - _START_TIME)

        # We're done, exit gracefully:
        sys.exit(0)
    except SystemExit as e:  # sys.exit()
        raise e
    except Exception as e:
        LOGGER.critical('ERROR: Unexpected Exception: %s', str(e))
        LOGGER.critical(str(e))
        traceback.print_exc()
        sys.exit(1)
