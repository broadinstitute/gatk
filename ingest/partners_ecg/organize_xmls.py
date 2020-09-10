import os
import re
import shutil
import logging
import datetime
import argparse
import xmltodict
from timeit import default_timer as timer


def _format_date(input_date, day_flag, sep_char='-'):
    '''Format input date from `mm-dd-yyyy` into `yyyy-mm`, or `yyyy-mm-dd`
    This is the ISO 8601 standard for date'''
    date_iso = input_date[6:10] + sep_char + input_date[0:2]
    if day_flag:
        date_iso = date_iso + sep_char + input_date[3:5]
    return date_iso


def _process_args(args):
    if not os.path.exists(args.source_xml_folder):
        raise NameError(f"{args.source_xml_folder} does not exist")
    if not os.path.exists(args.destination_xml_folder):
        os.mkdir(args.destination_xml_folder)
    if not os.path.exists(args.bad_xml_folder):
        os.mkdir(args.bad_xml_folder)

    now_str = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')

    # This is temporary becuase I can't figure out how to import load_config from ml4h.logger
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s - %(module)s:%(lineno)d - %(levelname)s - %(message)s",
                        handlers=[
                            logging.FileHandler(f"{now_str}_organize_xmls_log.txt"),
                            logging.StreamHandler()
                            ])


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--source_xml_folder",
                        default="/data/partners_ecg/xml",
                        help="Path to directory containing source XMLs")

    parser.add_argument("--destination_xml_folder",
                        default="/data/partners_ecg/dst",
                        help="Path to dir to organize XMLs in yyyy-mm dirs")

    parser.add_argument("--bad_xml_folder",
                        default="/data/partners_ecg/xml_bad",
                        help="Path to directory in which to store malformed XMLs")

    parser.add_argument("--method", default="copy", choices=["copy", "move"], help="Copy or move files from src to dst/yyyy-mm. Default: copy")

    parser.add_argument("--verbose", action="store_true", help="Print more information as each file is processed. Default: False")

    args = parser.parse_args()
    _process_args(args)
    return args


def run(args):
    start_time = timer()

    logging.info(f"Source XML location: {args.source_xml_folder}")
    logging.info(f"Destination XML location: {args.destination_xml_folder}")
    logging.info(f"Bad XML location: {args.bad_xml_folder}")

    num_bad_encodings = 0
    num_parsing_err = 0
    num_processed = 0
    num_files = 0
    for root, dirs, files in os.walk(args.source_xml_folder):
        for filename in files:
            if filename.endswith(".xml"):
                num_files += 1

    # Define accepted XML encodings
    valid_encodings = {"UTF-8", "UTF-16", "ISO-8859-1"}

    # Loop through XML files in the directory
    for root, dirs, files in os.walk(args.source_xml_folder):
        for filename in files:

            # Skip non-XML files
            if not filename.endswith(".xml"):
                continue

            fpath = os.path.join(root, filename)

            logging.debug(f"Parsing {fpath}")

            # Declare logical flag to indicate valid XML
            valid_xml = True

            # Open XML as a Python Dictionary
            with open(fpath) as fd:

                # Parse first line of XML to assess encoding;
                # if invalid, replace with valid encoding
                xml_as_string = fd.read()

                # Regex to find "encoding="STUFF-HERE"" within
                # <?xml version="1.0" encoding="STUFF-HERE"?>""
                # by returning everything after "encoding="" and before ""
                # re.search(r"(?<=encoding\=\")(.*?)(?=\"\?)", xml_as_string).group(0)

                # 1. Look behind positive (?<=B)A finds A preceded by B
                #    Here, our encoding value is left-bound by "encoding=""
                # 2. Return one or more characters via reluctant (lazy) match
                # 3. Look ahead postive A(?=B) finds expression A where B follows.
                #    Here, our encoding value is right-bound by ""?>"

                verbosePattern = re.compile("""
                    (?<=encoding\=\")
                    (.*?)
                    (?=\"\?\>)
                    """, re.VERBOSE)

                # Extract XML encoding from first line of imported XML
                xml_encoding = re.search(verbosePattern, xml_as_string)

                if xml_encoding is None:
                    valid_xml = False
                else:
                    xml_encoding = xml_encoding.group(0)

                # If xml_encoding is not among the accepted XML encodings, fix it
                if xml_encoding not in valid_encodings or not valid_xml:
                    logging.debug(f"Bad XML encoding found: {xml_encoding}. Replacing with ISO-8859-1.")

                    # Replace the bad encoding in xml_as_string with ISO-8859-1
                    xml_as_string = re.sub(verbosePattern, "ISO-8859-1", xml_as_string, count=1)

                    # Increment counter to track bad encodings
                    num_bad_encodings += 1

                    # Overwrite XML file with corrected encoding
                    with open(fpath, "w") as f:
                        f.write(xml_as_string)

                try:
                    # Parse XMl-as-string into a dict
                    doc = xmltodict.parse(xml_as_string)

                    # Isolate acquisition date of test and save as mm-dd-yyyy format
                    ecg_date = doc["RestingECG"]["TestDemographics"]["AcquisitionDate"]

                    logging.debug("Acquisition date found! " + ecg_date)

                    # Check if the date
                    # 1) has ten digits
                    # 2) has a dash at index 2
                    # 3) has a dash at index 5
                    date_check = [len(ecg_date) == 10,
                                  ecg_date[2] == "-",
                                  ecg_date[5] == "-"]

                    # If does not pass date check, inform user XML has bad date format
                    if not all(date_check):
                        valid_xml = False

                    # If passes date check
                    else:

                        # Define full path to new directory in yyyy-mm format
                        yyyymm = _format_date(ecg_date, day_flag=False)
                        args.dst_yyyymm = os.path.join(args.destination_xml_folder, yyyymm)

                        # If directory does not exist, create it
                        if not os.path.exists(args.dst_yyyymm):
                            os.makedirs(args.dst_yyyymm)
                            logging.debug(f"No valid yyyy-mm directory exists. Creating: {args.dst_yyyymm}")

                # If there is any parsing error, set flag to False
                except xmltodict.expat.ExpatError:
                    valid_xml = False
                    num_parsing_err += 1

            # If the XML is not valid, set the final path to bad xml path
            if not valid_xml:
                args.dst_yyyymm = args.bad_xml_folder
                logging.debug("Missing or invalid acquisition date, or other XML error.")

            # If new directory does not exist, create it
            if not os.path.exists(args.dst_yyyymm):
                os.makedirs(args.dst_yyyymm)

            # Copy or move XML file into new directory
            fpath_xml_newdir = args.dst_yyyymm + "/" + filename
            if 'copy' == args.method:
                shutil.copy(fpath, fpath_xml_newdir)
                logging.debug(f"XML copied to {args.dst_yyyymm}")
            elif 'move' == args.method:
                shutil.move(fpath, fpath_xml_newdir)
                logging.debug(f"XML moved to {args.dst_yyyymm}")

            num_processed += 1

            # Log progress every 100 files
            if num_processed % 100 == 0:
                logging.info(f"processed {num_processed} / {num_files} XML files ({num_processed/num_files*100:.1f}% done)")

    # Log final results
    num_valid = num_processed - num_bad_encodings - num_parsing_err
    logging.info(f"Number files found in src: {num_files}")
    logging.info(f"Number valid files now in dst: {num_valid}")
    logging.info(f"Number files with bad encodings or parsing errors: {num_bad_encodings + num_parsing_err}")

    end_time = timer()
    elapsed_time = end_time - start_time
    logging.info(f"Executed organize script in {elapsed_time:.2f} sec")


if __name__ == "__main__":
    args = parse_args()
    run(args)
