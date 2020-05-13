import os
import re
import array
import base64
import struct
import logging
from datetime import datetime
from abc import ABC, abstractmethod
from collections import defaultdict
from typing import List, Dict, Tuple, Union
from xml.parsers.expat import ExpatError, ParserCreate

import h5py
import numcodecs
import numpy as np
from joblib import Parallel, delayed
from bs4 import BeautifulSoup, SoupStrainer, PageElement

from ml4cvd.defines import TENSOR_EXT, XML_EXT


ECG_REST_INDEPENDENT_LEADS = ["I", "II", "V1", "V2", "V3", "V4", "V5", "V6"]


def write_tensors_partners(xml_folder: str, tensors: str) -> None:
    """Write tensors as HD5 files containing data from Partners dataset

    One HD5 is generated per patient. One HD5 may contain multiple ECGs.

    :param xml_folder: Path to folder containing ECG XML files organized in
                       subfolders by date
    :param tensors: Folder to populate with HD5 tensors

    :return: None
    """

    n_jobs = -1

    logging.info('Mapping XMLs to MRNs')
    mrn_xmls_map = _get_mrn_xmls_map(xml_folder, n_jobs=n_jobs)

    logging.info('Converting XMLs into HD5s')
    _convert_mrn_xmls_to_hd5_wrapper(mrn_xmls_map, tensors, n_jobs=n_jobs)


def _map_mrn_to_xml(fpath_xml: str) -> Union[Tuple[str, str], None]:
    with open(fpath_xml, 'r') as f:
        for line in f:
            match = re.match(r'.*<PatientID>(.*)</PatientID>.*', line)
            if match:
                mrn = _clean_mrn(match.group(1), fallback='bad_mrn')
                return (mrn, fpath_xml)
    logging.warning(f'No PatientID found at {fpath_xml}')
    return None


def _get_mrn_xmls_map(xml_folder: str, n_jobs: int) -> Dict[str, List[str]]:

    # Get all xml paths
    fpath_xmls = []
    for root, dirs, files in os.walk(xml_folder):
        for file in files:
            if os.path.splitext(file)[-1].lower() != XML_EXT:
                continue
            fpath_xmls.append(os.path.join(root, file))
    logging.info(f'Found {len(fpath_xmls)} XMLs at {xml_folder}')

    # Read through xmls to get MRN in parallel
    mrn_xml_list = Parallel(n_jobs=n_jobs)(delayed(_map_mrn_to_xml)(fpath_xml) for fpath_xml in fpath_xmls)

    # Build dict of MRN to XML files with that MRN
    mrn_xml_dict = defaultdict(list)
    for mrn_xml in mrn_xml_list:
        if mrn_xml:
            mrn_xml_dict[mrn_xml[0]].append(mrn_xml[1])
    logging.info(f'Found {len(mrn_xml_dict)} distinct MRNs')

    return mrn_xml_dict


def _clean_mrn(mrn: str, fallback: str) -> str:
    # TODO additional cleaning like o->0, |->1
    try:
        clean = re.sub(r'[^0-9]', '', mrn)
        clean = int(clean)
        if not clean:
            raise ValueError()
        return str(clean)
    except ValueError:
        logging.warning(f'Could not clean MRN "{mrn}" to an int. Falling back to "{fallback}".')
        return fallback


def _clean_read_text(text: str) -> str:
    # Convert to lowercase
    text = text.lower()

    # Replace newline character with space
    text = re.sub(r'\n', ' ', text)

    # Remove punctuation
    text = re.sub(r'[^\w\s]', '', text)

    # Replace two+ spaces with one space
    text = re.sub(r'  +', ' ', text)

    # Remove all leading and trailing whitespace
    text = text.strip()

    return text


def _get_voltage_from_xml(fpath_xml: str) -> Tuple[Dict[str, np.ndarray], str]:
    g_parser = MuseXmlParser()

    def start_element(name, attrs):
        g_parser.start_element(name, attrs)

    def end_element(name):
        g_parser.end_element(name)

    def char_data(data):
        g_parser.char_data(data)

    p = ParserCreate()
    p.StartElementHandler = start_element
    p.EndElementHandler = end_element
    p.CharacterDataHandler = char_data

    # Read the XML file and parse it
    p.ParseFile(open(fpath_xml, 'rb'))

    # Convert the data into a ZCG buffer
    g_parser.makeZcg()

    # Return voltage data as a dict of arrays
    voltage, units = g_parser.get_voltage()

    # Return voltage
    return voltage, units


def _text_from_xml(fpath_xml: str) -> Dict[str, str]:
    # Initialize empty dictionary in which to store text from the XML
    ecg_data = dict()

    # Define tags that we want to find and use SoupStrainer to speed up search
    tags = [
        'patientdemographics',
        'testdemographics',
        'order',
        'restingecgmeasurements',
        'originalrestingecgmeasurements',
        'intervalmeasurementtimeresolution',
        'intervalmeasurementamplituderesolution',
        'intervalmeasurementfilter',
        'diagnosis',
        'originaldiagnosis',
    ]
    strainer = SoupStrainer(tags)

    # Use lxml parser, which makes all tags lower case
    with open(fpath_xml, 'r') as f:
        soup = BeautifulSoup(f, 'lxml', parse_only=strainer)

    # If the XML is gibberish and un-parseable,
    # then soup.prettify() will return '' or False
    # and we should return an empty dict()

    # If the XML is not gibberish, parse the contents
    if soup.prettify():

        # Loop through the tags we want to extract
        for tag in tags:
            append_tag = ''
            if tag == 'restingecgmeasurements':
                append_tag = '_md'
            elif tag == 'originalrestingecgmeasurements':
                append_tag = '_pc'
            elif tag == 'diagnosis':
                # Parse text of cardiologist read within <Diagnosis> tag and save to ecg_data dict
                ecg_data['diagnosis_md'] = _parse_soup_diagnosis(soup.find('diagnosis'))
                continue
            elif tag == 'originaldiagnosis':
                # Parse text of computer read within <OriginalDiagnosis> tag and save to ecg_data dict
                ecg_data['diagnosis_pc'] = _parse_soup_diagnosis(soup.find('originaldiagnosis'))
                continue

            # Ovewrite soup with contents in major tag
            soup_of_tag = soup.find(tag)

            # Check if soup_of_tag is not None, otherwise find_all() throws an
            # error if called on a 'NoneType' object
            if soup_of_tag is not None:
                # Find all child elements in subset
                elements = soup_of_tag.find_all()

                # If there are no child elements, get the text of the tag
                if not elements:
                    elements = [soup_of_tag]

                # Iterate through child elements via list comprehension
                # and save the element key-value (name-text) to dict
                [ecg_data.update({el.name + append_tag: el.get_text()}) for el in elements]

    # Return dict
    return ecg_data


def _parse_soup_diagnosis(input_from_soup: PageElement) -> str:

    parsed_text = ''

    if input_from_soup is None:
        return parsed_text
    else:
        parts = input_from_soup.find_all('diagnosisstatement')

        # Check for edge case where <diagnosis> </diagnosis> does not encompass
        # <DiagnosisStatement> sub-element, which results in parts being length 0
        if len(parts) > 0:
            for part in parts:

                # Create list of all <stmtflag> entries
                flags = part.find_all('stmtflag')

                # Isolate text from part
                text_to_append = part.find('stmttext').text

                # Initialize flag to ignore sentence, e.g. do not append it
                flag_ignore_sentence = False

                # If no reasons found, append
                if not flag_ignore_sentence:
                    # Append diagnosis string with contents within <stmttext> tags
                    parsed_text += text_to_append

                    endline_flag = False

                    # Loop through flags and if 'ENDSLINE' found anywhere, mark flag
                    for flag in flags:
                        if flag.text == 'ENDSLINE':
                            endline_flag = True

                    # If 'ENDSLINE' was found anywhere, append newline
                    if endline_flag:
                        parsed_text += '\n'

                    # Else append space
                    else:
                        parsed_text += ' '

            # Remove final newline character in diagnosis
            if parsed_text[-1] == '\n':
                parsed_text = parsed_text[:-1]

        return parsed_text


def _compress_and_save_data(
    hd5: h5py.Group, name: str, data: Union[str, np.ndarray],
    dtype: str, method: str = 'zstd', compression_opts: int = 19,
) -> None:
    # Define codec
    codec = numcodecs.zstd.Zstd(level=compression_opts)

    # If data is string, encode to bytes
    if dtype == 'str':
        data_compressed = codec.encode(data.encode())
        dsize = len(data.encode())
    else:
        data_compressed = codec.encode(data)
        dsize = len(data) * data.itemsize

    # Save data to HD5
    dat = hd5.create_dataset(name=name, data=np.void(data_compressed))

    # Set attributes
    dat.attrs['method'] = method
    dat.attrs['compression_level'] = compression_opts
    dat.attrs['len'] = len(data)
    dat.attrs['uncompressed_length'] = dsize
    dat.attrs['compressed_length'] = len(data_compressed)
    dat.attrs['dtype'] = dtype


def _get_max_voltage(voltage: Dict[str, np.ndarray]) -> float:
    max_voltage = 0
    for lead in voltage:
        if max(voltage[lead]) > max_voltage:
            max_voltage = max(voltage[lead])
    return max_voltage


def _convert_xml_to_hd5(fpath_xml: str, fpath_hd5: str, hd5: h5py.Group) -> int:
    # Return 1 if converted, 0 if ecg was bad or -1 if ecg was a duplicate
    # Set flag to check if we should convert to hd5
    convert = 1

    # Extract text data from XML into dict
    text_data = _text_from_xml(fpath_xml)
    dt = datetime.strptime(f"{text_data['acquisitiondate']} {text_data['acquisitiontime']}", '%m-%d-%Y %H:%M:%S')
    ecg_dt = dt.isoformat()

    # If XML is empty, remove the XML file and set convert to false
    if (os.stat(fpath_xml).st_size == 0 or not text_data):
        os.remove(fpath_xml)
        convert = 0
        logging.warning(f'Conversion of {fpath_xml} failed! XML is empty.')
    else:
        # Check if patient already has an ECG at given date and time
        if ecg_dt in hd5.keys():
            logging.warning(f'Conversion of {fpath_xml} skipped. Converted XML already exists in HD5.')
            convert = -1

        # Extract voltage from XML
        try:
            voltage, voltage_units = _get_voltage_from_xml(fpath_xml)

            # If the max voltage value is 0, do not convert
            if _get_max_voltage(voltage) == 0:
                logging.warning(f'Conversion of {fpath_xml} failed! Maximum voltage is 0.')
                convert = 0

        # If there is no voltage, or the XML is poorly formed,
        # the function will throw an exception, and we mark 'convert' to False
        # However, ExpatError should be impossible to throw, since earlier
        # we catch XMLs filled with gibberish
        except (IndexError, ExpatError):
            logging.warning(f'Conversion of {fpath_xml} failed! Voltage is empty or badly formatted.')
            convert = 0

    # If convert is still true up to here, make the hd5
    if convert == 1:

        # Define keys for cleaned reads
        key_read_md = 'diagnosis_md'
        key_read_pc = 'diagnosis_pc'
        key_read_md_clean = 'read_md_clean'
        key_read_pc_clean = 'read_pc_clean'

        # Clean cardiologist read
        read_md_clean = None
        if key_read_md in text_data.keys():
            read_md_clean = _clean_read_text(text=text_data[key_read_md])

        # Clean MUSE read
        read_pc_clean = None
        if key_read_pc in text_data.keys():
            read_pc_clean = _clean_read_text(text=text_data[key_read_pc])

        # Clean Patient MRN to only numeric characters
        key_mrn_clean = 'patientid_clean'
        mrn_clean = None
        if 'patientid' in text_data.keys():
            mrn_clean = _clean_mrn(text_data['patientid'], fallback='')

        gp = hd5.create_group(ecg_dt)

        # Save cleaned MRN to hd5
        if mrn_clean:
            _compress_and_save_data(hd5=gp, name=key_mrn_clean, data=mrn_clean, dtype='str')

        # Iterate through each lead and save voltage array to hd5
        for lead in voltage.keys():
            _compress_and_save_data(hd5=gp, name=lead, data=voltage[lead].astype('int16'), dtype='int16')

        # Save voltage units to hd5
        key_voltage_units = 'voltageunits'
        if voltage_units != '':
            _compress_and_save_data(hd5=gp, name=key_voltage_units, data=voltage_units, dtype='str')

        # Iterate through keys in extracted text data and save to hd5
        for key in text_data.keys():
            _compress_and_save_data(hd5=gp, name=key, data=text_data[key], dtype='str')

        # Save cleaned reads to hd5
        if read_md_clean:
            _compress_and_save_data(hd5=gp, name=key_read_md_clean, data=read_md_clean, dtype='str')

        if read_pc_clean:
            _compress_and_save_data(hd5=gp, name=key_read_pc_clean, data=read_pc_clean, dtype='str')

        logging.info(f'Wrote {fpath_xml} to {fpath_hd5}')
    return convert


def _convert_mrn_xmls_to_hd5(mrn: str, fpath_xmls: List[str], dir_hd5: str, hd5_prefix: str) -> Tuple[int, int, int]:
    fpath_hd5 = os.path.join(dir_hd5, f'{mrn}{TENSOR_EXT}')
    num_xml_converted = 0
    num_dupe_skipped = 0
    num_src_in_hd5 = 0
    num_ecg_in_hd5 = 0

    with h5py.File(fpath_hd5, 'a') as hd5:
        hd5_ecg = hd5[hd5_prefix] if hd5_prefix in hd5.keys() else hd5.create_group(hd5_prefix)
        for fpath_xml in fpath_xmls:
            converted = _convert_xml_to_hd5(fpath_xml, fpath_hd5, hd5_ecg)
            if converted == 1:
                num_xml_converted += 1
            elif converted == -1:
                num_dupe_skipped += 1
        num_ecg_in_hd5 = len(hd5_ecg.keys())

        # If there are no ECGs in HD5, delete ECG group
        # There may be prior ECGs in HD5
        # num_xml_converted != num_ecg_in_hd5
        if not num_ecg_in_hd5: del hd5[hd5_prefix]

        num_src_in_hd5 = len(hd5.keys())

    if not num_src_in_hd5:
        # If there is no other data in HD5, delete HD5
        try:
            os.remove(fpath_hd5)
        except:
            logging.warning(f'Could not delete empty HD5 at {fpath_hd5}')

    num_hd5_written = 1 if num_xml_converted else 0

    return (num_hd5_written, num_xml_converted, num_dupe_skipped)


def _convert_mrn_xmls_to_hd5_wrapper(mrn_xmls_map: Dict[str, List[str]], dir_hd5: str, hd5_prefix: str = 'partners_ecg_rest', n_jobs: int = -1):
    tot_xml = sum([len(v) for k, v in mrn_xmls_map.items()])
    os.makedirs(dir_hd5, exist_ok=True)
    converted = Parallel(n_jobs=n_jobs)(delayed(_convert_mrn_xmls_to_hd5)(mrn, fpath_xmls, dir_hd5, hd5_prefix) for mrn, fpath_xmls in mrn_xmls_map.items())
    num_hd5 = sum([x[0] for x in converted])
    num_xml = sum([x[1] for x in converted])
    num_dup = sum([x[2] for x in converted])

    logging.info(f"Converted {num_xml} XMLs to {num_hd5} HD5s at {dir_hd5}")
    logging.info(f"Skipped {num_dup} duplicate XMLs")
    logging.info(f"Skipped {tot_xml - num_dup - num_xml} malformed XMLs")


class XmlElementParser(ABC):
    """Abstract base class for a XML Parsing State. It contains methods for
    restoring the previous state and for tracking the character data between
    tags"""
    def __init__(self, old_State=None):
        self.__old_State = old_State
        self.__data_Text = ""
        super().__init__()

    def restoreState(self, context):
        """This method restores the previous state in the XML parser"""
        if self.__old_State:
            context.setState(self.__old_State)

    def clearData(self):
        """This method clears the character data that was collected during parsing"""
        self.__data_Text = ""

    def getData(self):
        """This method returns the character data that was collected during
        parsing and it strips any leading or trailing whitespace"""
        return str.strip(self.__data_Text)

    @abstractmethod
    def start_element(self, name, attrs, context):
        pass

    @abstractmethod
    def end_element(self, name, context):
        pass

    def char_data(self, data, context):
        """This method accumulates any character data"""
        self.__data_Text = self.__data_Text + data


class IdleParser(XmlElementParser):
    """State for handling the Idle condition"""
    def __init__(self):
        XmlElementParser.__init__(self)

    def start_element(self, name, attrs, context):
        if name == WaveformElementParser.Tag:
            context.setState(WaveformElementParser(self))

    def end_element(self, name, context):
        self.clearData()


class WaveformElementParser(XmlElementParser):
    """State for handling the Waveform element"""

    Tag = "Waveform"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()
        if name in (ElementParser.waveform_type_tag, ElementParser.voltage_base_tag):
            context.setState(ElementParser(self))
        elif name == LeadDataElementParser.Tag:
            context.setState(LeadDataElementParser(self))

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)


class LeadDataElementParser(XmlElementParser):
    """State for handling the LeadData element"""

    Tag = "LeadData"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()
        if name in (
            ElementParser.lead_id_tag, ElementParser.unit_tag,
            ElementParser.bit_tag, ElementParser.waveform_data_tag,
        ):
            context.setState(ElementParser(self))

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)


class ElementParser(XmlElementParser):
    """State for handling voltageBase Element, LeadID,
    LeadAmplitudeUnits, LeadAmplitudeUnitsPerBit, WaveformData
    and WaveFormType elements"""

    voltage_base_tag = "voltageBase"
    unit_tag = "LeadAmplitudeUnits"
    bit_tag = "LeadAmplitudeUnitsPerBit"
    lead_id_tag = "LeadID"
    waveform_type_tag = "WaveformType"
    waveform_data_tag = "WaveFormData"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.voltage_base_tag:
            self.restoreState(context)
            if context.found_Rhythm:
                context.setvoltageBase(self.getData())
        if name == self.unit_tag:
            self.restoreState(context)
            context.setUnits(str.strip(self.getData()))
        if name == self.bit_tag:
            self.restoreState(context)
            context.setAdu(float(str.strip(self.getData())))
        if name == self.lead_id_tag:
            self.restoreState(context)
            if context.found_Rhythm:
                context.addLeadId(self.getData())
        if name == self.waveform_type_tag:
            self.restoreState(context)
            if str.find(self.getData(), "Rhythm") >= 0:
                context.setRhythmFound(1)
            else:
                context.setRhythmFound(0)
        if name == self.waveform_data_tag:
            self.restoreState(context)
            if context.found_Rhythm:
                context.addWaveformData(self.getData())


class MuseXmlParser:
    """This class is the parsing context in the object-oriented State pattern"""
    def __init__(self):
        self.ecg_Data = dict()
        self.ecg_leads = list()
        self.__state = IdleParser()
        self.found_Rhythm = 0
        self.voltage_Rate = 0
        self.adu_Gain = 1
        self.units = ""

    def setState(self, s):
        self.__state = s

    def getState(self):
        return self.__state

    def setvoltageBase(self, text):
        if self.found_Rhythm:
            self.voltage_Rate = int(text)

    def setAdu(self, gain):
        self.adu_Gain = gain

    def setUnits(self, units):
        self.units = units

    def setRhythmFound(self, v):
        self.found_Rhythm = v

    def addLeadId(self, id):
        self.lead_Id = id

    def addWaveformData(self, text):
        self.ecg_Data[self.lead_Id] = base64.b64decode(text)
        self.ecg_leads.append(self.lead_Id)

    def start_element(self, name, attrs):
        """This function trackes the start elements found in the XML file with a
        simple state machine"""
        self.__state.start_element(name, attrs, self)

    def end_element(self, name):
        self.__state.end_element(name, self)

    def char_data(self, data):
        self.__state.char_data(data, self)

    def makeZcg(self):
        """This function converts the data read from the XML file into a ZCG buffer
        suitable for storage in binary format"""

        # All of the leads should have the same number of voltage
        total_byte_chars = len(self.ecg_Data[self.ecg_leads[0]])

        # We have 2 bytes per ECG voltage, so make our buffer size n * DATAMUX
        self.zcg = array.array('d')

        # Append the data into our huge ZCG buffer in the correct order
        for t in range(0, total_byte_chars, 2):
            for lead in self.ecg_leads:
                voltage = struct.unpack("h", bytes([self.ecg_Data[lead][t], self.ecg_Data[lead][t+1]]))
                self.zcg.append(voltage[0])

    def get_voltage(self):
        # Multiply base voltage by the gain
        self.zcg = [x * self.adu_Gain for x in self.zcg]

        # Initialize dict of lists where L keys are ECG leads
        # L is not always 8; sometimes you have right-sided leads (e.g. V3R)
        voltage = {new_list: [] for new_list in self.ecg_leads}

        # Loop through all indices in 1D voltage array, L leads at a time
        # Correctly parsing the 1D array requires we parse all L leads
        for i in range(0, len(self.zcg), len(self.ecg_leads)):
            # Initialize lead counter; this is necessary to increment through
            # the 1D representation of a 2D array. For example, we start at
            # element i=0, then loop through element i+1, i+2, etc. until we
            # iterate through all 8 leads. Then we move to element i=8 (which
            # is the step size in the outermost for loop), and repeat.
            #
            #  0 [5 ... ]
            #  1 [7 ... ]
            #  2 [2 ... ]
            #  3 [3 ... ]
            #  4 [1 ... ]
            #  5 [1 ... ]
            #  6 [4 ... ]
            #  7 [3 ... ]

            k = 0

            for lead in self.ecg_leads:

                voltage[lead].append(self.zcg[i + k])
                k = k + 1

        # Overwrite voltage dict with only standard 8 leads
        voltage = {lead: voltage[lead] for lead in ECG_REST_INDEPENDENT_LEADS}

        # Convert dict of lists into dict of arrays to enable vector math
        voltage = {key: np.array(val) for key, val in voltage.items()}

        # Perform vector math to obtain remaining 4 leads
        voltage['III'] = voltage['II'] - voltage['I']
        voltage['aVR'] = -1 * (voltage['I'] + voltage['II']) / 2
        voltage['aVL'] = voltage['I'] - voltage['II'] / 2
        voltage['aVF'] = voltage['II'] - voltage['I'] / 2

        return voltage, self.units
