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
import itertools
from collections import namedtuple
from tqdm import tqdm

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

CONTIG_PCT_POS_CHECK = 30

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

    _parser.add_argument('-d', "--detailed", action='store_true',
                         help="Perform detailed all-to-all base-level analysis on all contigs that are not identical by"
                              " sequence dictionary md5sum."
                              "This will check lengths of all contigs.  Then for all equal length contigs a portion of"
                              " the contig will be examined approxmiately 30% of the way into the contig.  If these "
                              "subsequences are identical or nearly identical then the contigs will be processed and "
                              "exact differences will be emitted as variants.")

    _parser.add_argument('--case_sensitive', action='store_true',
                         help="Set the comparisons between bases between references to be case-sensitive.  "
                              "(default: case insensitive)")

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
ContigBaseInfo = namedtuple("ContigBaseInfo", ['ref', 'contig', 'pos', 'base'])


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


def get_ref_base_name(ref_fasta_path):
    """Get the base name to be used to represent the given reference FASTA file."""
    base_path, _ = os.path.splitext(ref_fasta_path)
    return os.path.basename(base_path)


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
        seq_dict = get_sequence_dict(base_path + ".dict")

        base_name = get_ref_base_name(fasta_path)
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
    printed_md5s = set()
    for ref_name in sorted_ref_names:
        identical_contig_list = sorted(identical_contig_map[ref_name].keys(), key=natural_keys)
        for contig_name in identical_contig_list:
            contig_md5 = seq_dict_map[ref_name][contig_name].md5sum
            if contig_md5 in printed_md5s:
                continue

            fields_to_print = [contig_md5]
            for ref in sorted_ref_names:
                if ref == ref_name:
                    fields_to_print.append(contig_name)
                elif ref in identical_contig_map[ref_name][contig_name]:
                    fields_to_print.append(identical_contig_map[ref_name][contig_name][ref])
                else:
                    fields_to_print.append("----")

            print(row_format_string.format(*fields_to_print))
            printed_md5s.add(contig_md5)

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


def sample_and_check_contig_seqs(ref_files, ref_contig_dict, contig_length, num_bases_to_sample=50):
    """Samples a number of bases on each given reference / contig pair.
    Bases are sampled at +/- CONTIG_PCT_POS_CHECK % of the length of the contig.
    The comparison here is always case INSENSITIVE.
    If either the "Front" or "Back" bases are the same between all given ref/contig pairs, then this function returns
    True.  Otherwise returns false."""

    front_base_samples = set()
    back_base_samples = set()

    for fasta_path in ref_files:
        base_name = get_ref_base_name(fasta_path)
        sample_start = int((CONTIG_PCT_POS_CHECK/100.0) * contig_length)
        sample_start_back = contig_length - sample_start

        contig = ref_contig_dict[base_name]

        with pysam.FastaFile(fasta_path) as f:
            front_base_samples.add(f.fetch(contig, sample_start, sample_start + num_bases_to_sample).upper())
            back_base_samples.add(f.fetch(contig, sample_start_back, sample_start_back + num_bases_to_sample).upper())

    return len(front_base_samples) == 1 or len(back_base_samples) == 1


def scrutinize_contig(ref_files, ref_contig_dict, contig_length, large_diff_thresh=None, ignore_case=True):
    """Iterates over all bases in the given contigs tabulating differences between them.
    If large differences are discovered, will stop comparison."""

    differences = dict()
    open_fasta_files = []

    try:
        num_differences = 0

        base_names = []
        contig_names = []
        for fasta_path in ref_files:
            open_fasta_files.append(pysam.FastaFile(fasta_path))

            base_name = get_ref_base_name(fasta_path)
            base_names.append(base_name)
            contig_names.append(ref_contig_dict[base_name])

        step_size = 1024

        start = 0
        end = start + step_size
        if end > contig_length - 1:
            end = contig_length

        # NOTE: Interval conventions for pysam are 0-based, half-open
        with tqdm(total=contig_length, desc="Contig: " + contig_names[0], unit=" bases") as pbar:
            while True:
                bases_list = []
                bases_hash_set = set()
                for i, f in enumerate(open_fasta_files):
                    if ignore_case:
                        b = f.fetch(contig_names[i], start, end).upper()
                    else:
                        b = f.fetch(contig_names[i], start, end)
                    bases_list.append(b)
                    bases_hash_set.add(hash(b))

                # Do we have any differences in our large chunk?
                if len(bases_hash_set) > 1:
                    for base_num in range(len(bases_list[0])):
                        # Get the bases that could be different:
                        b_set = set()
                        b_description_list = []
                        pos = base_num + start
                        for ref_num in range(len(bases_list)):

                            if ignore_case:
                                b = bases_list[ref_num][base_num].upper()
                            else:
                                b = bases_list[ref_num][base_num]

                            b_set.add(b)
                            b_description_list.append(
                                ContigBaseInfo(base_names[ref_num],
                                               contig_names[ref_num],
                                               pos,
                                               bases_list[ref_num][base_num])
                            )
                        if len(b_set) != 1:
                            differences[pos] = b_description_list
                            num_differences += 1

                        # Update our progress bar by 1 base:
                        pbar.update(1)

                        if large_diff_thresh and num_differences > large_diff_thresh:
                            LOGGER.warning("Found a lot of differences.  Halting comparison.")
                            return differences
                else:
                    # Update our progress bar by our chunk size:
                    pbar.update(step_size)

                # Update our position and exit if necessary:
                if end == contig_length:
                    break

                start = end
                end = start + step_size
                if end > contig_length - 1:
                    end = contig_length
    finally:
        for f in open_fasta_files:
            if not f.closed:
                f.close()

    return differences


def get_scrutinized_row_format_string(differences, ref_names, col_spacer='\t'):
    """Get the format string for each row in the given scrutinized differences."""

    # Output should be like this:
    # REF_1.contig_name REF_1.pos REF_1.base REF_2.contig_name REF_2.pos REF_2.base [REF_3.contig_name ...]

    ref_col_max_sizes = dict()
    for ref in ref_names:
        ref_col_max_sizes[ref] = {'contig': 6, "pos": 3}

    for contig_diffs in differences:
        for ref_diff in contig_diffs.values():
            for row in ref_diff:
                ref = row.ref
                if len(row.contig) > ref_col_max_sizes[ref]['contig']:
                    ref_col_max_sizes[ref]['contig'] = len(row.contig)
                if len(str(row.pos)) > ref_col_max_sizes[ref]['pos']:
                    ref_col_max_sizes[ref]['pos'] = len(str(row.pos))

    row_format_string = ""
    for ref in ref_names:
        row_format_string = row_format_string + \
                            "{:<" + str(ref_col_max_sizes[ref]['contig'] + len(ref) + 1) + "}" + col_spacer + \
                            "{:<" + str(ref_col_max_sizes[ref]['pos'] + len(ref) + 1) + "}" + col_spacer + \
                            "{:<" + str(4 + len(ref) + 1) + "}" + col_spacer

    # Remove the last spacer:
    row_format_string = row_format_string[:-len(col_spacer)]

    return row_format_string


def print_scrutinized_table_header(differences, ref_names, col_spacer='\t'):
    """Prints the header for a master scrutinized contig table."""
    format_string = get_scrutinized_row_format_string(differences, ref_names, col_spacer)

    row = list(itertools.chain(*[[f"{ref}.contig", f"{ref}.pos", f"{ref}.base"] for ref in ref_names]))
    print(format_string.format(*row))


def print_scrutinized_table(differences, ref_names, col_spacer='\t'):
    """Prints a table of the scrutinized contig diff information."""
    print_scrutinized_table_header(differences, ref_names, col_spacer)

    format_string = get_scrutinized_row_format_string(differences, ref_names, col_spacer)

    for contig_diffs in differences:
        for ref_diff in contig_diffs.values():
            ordered_row_fields = []
            for ref in ref_names:
                for row in ref_diff:
                    if ref != row.ref:
                        continue
                    ordered_row_fields.append(row.contig)
                    ordered_row_fields.append(row.pos)
                    ordered_row_fields.append(row.base)
            print(format_string.format(*ordered_row_fields))


def perform_detailed_analysis(_args, identical_contig_map, seq_dict_map):
    """Perform detailed all-to-all base-level analysis on all contigs that are not identical by
    sequence dictionary md5sum.

    This will check lengths of all contigs.  Then for all equal length contigs a portion of
    the contig will be examined approxmiately 30% of the way into the contig.  If these
    subsequences are identical or nearly identical then the contigs will be processed and
    exact differences will be emitted as variants."""

    num_refs = len(seq_dict_map.keys())
    LOGGER.info("Performing detailed analysis of %d references...", num_refs)

    contig_diff_lists = []

    # First get a list of contigs that might be equal:
    contig_lengths_dict = dict()
    for ref in seq_dict_map.keys():
        for contig in seq_dict_map[ref]:
            contig_length = seq_dict_map[ref][contig].length
            if contig_length not in contig_lengths_dict:
                contig_lengths_dict[contig_length] = {ref: contig}
            else:
                contig_lengths_dict[contig_length][ref] = contig

    # Remove any contigs that have unique lengths that are
    # represented across all references so we do less work:
    contig_lengths = list(contig_lengths_dict.keys())
    lengths_to_remove = set()
    for l in contig_lengths:
        if len(contig_lengths_dict[l]) == 1:
            lengths_to_remove.add(l)
            continue
        for ref, contig_name in contig_lengths_dict[l].items():
            if contig_name in identical_contig_map[ref]:
                if len(identical_contig_map[ref][contig_name]) == (num_refs - 1):
                    lengths_to_remove.add(l)

    for l in lengths_to_remove:
        del contig_lengths_dict[l]

    if LOGGER.isEnabledFor(logging.INFO):
        LOGGER.info("Possibly identical seqs: ")
        for l in contig_lengths_dict.keys():
            LOGGER.info("    Length: %d: %s", l,
                        ", ".join([f"{ref}[{contig_name}]" for ref, contig_name in contig_lengths_dict[l].items()]))

    LOGGER.info("Sampling bases from %d%% through contigs to determine if full check should be performed.",
                CONTIG_PCT_POS_CHECK)

    # Now process the remaining contigs:
    for l in contig_lengths_dict.keys():
        if not sample_and_check_contig_seqs(_args.ref_files, contig_lengths_dict[l], l):
            LOGGER.info("Failed filtering: %d - %s", l, [f"{k}[{v}]" for k, v in contig_lengths_dict[l].items()])
        else:
            LOGGER.info("Passed filtering:  %d - %s", l, [f"{k}[{v}]" for k, v in contig_lengths_dict[l].items()])
            differences = scrutinize_contig(_args.ref_files, contig_lengths_dict[l], l,
                                            ignore_case=not _args.case_sensitive)
            LOGGER.info(f"Found {len(differences)} differences on contig: "
                        f"{' / '.join(contig_lengths_dict[l].values())}")

            for ref_differences in differences.values():
                fields_to_print = [ref_differences[0].contig, str(ref_differences[0].pos)]
                for d in ref_differences:
                    fields_to_print.append(d.base)

            contig_diff_lists.append(differences)

    return contig_diff_lists

########################################
# Main function:
####################


def _main(_args):
    """Default main implementation."""

    # Do a dictionary comparison:
    unique_contig_map, identical_contig_map, seq_dict_map = analyze_sequence_dictionaries(_args)
    print_contig_table(unique_contig_map, identical_contig_map, seq_dict_map)

    # Do an in-depth comparison of the non-identical contigs if requested:
    if _args.detailed:
        print("=" * 80)
        print("=" * 80)
        print("=" * 80)

        differences = perform_detailed_analysis(_args, identical_contig_map, seq_dict_map)
        print_scrutinized_table(differences, [get_ref_base_name(f) for f in _args.ref_files])


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

