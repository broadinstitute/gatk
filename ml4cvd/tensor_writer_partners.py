import os
import re
import base64
import struct
import logging
import multiprocessing
from datetime import datetime
from collections import defaultdict
from typing import List, Dict, Tuple, Union

import bs4
import h5py
import numcodecs
import numpy as np

from ml4cvd.defines import TENSOR_EXT, XML_EXT, ECG_REST_AMP_LEADS

ECG_REST_INDEPENDENT_LEADS = ['I', 'II', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6']


def write_tensors_partners(xml_folder: str, tensors: str, num_workers: int) -> None:
    """Write tensors as HD5 files containing data from Partners dataset

    One HD5 is generated per patient. One HD5 may contain multiple ECGs.

    :param xml_folder: Path to folder containing ECG XML files organized in
                       subfolders by date
    :param tensors: Folder to populate with HD5 tensors

    :return: None
    """

    logging.info('Mapping XMLs to MRNs')
    mrn_xmls_map = _get_mrn_xmls_map(xml_folder, num_workers)

    logging.info('Converting XMLs into HD5s')
    _convert_mrn_xmls_to_hd5_wrapper(mrn_xmls_map, tensors, num_workers)


def _map_mrn_to_xml(fpath_xml: str) -> Union[Tuple[str, str], None]:
    with open(fpath_xml, 'r') as f:
        for line in f:
            match = re.match(r'.*<PatientID>(.*)</PatientID>.*', line)
            if match:
                mrn = _clean_mrn(match.group(1), fallback='bad_mrn')
                return (mrn, fpath_xml)
    logging.warning(f'No PatientID found at {fpath_xml}')
    return None


def _get_mrn_xmls_map(xml_folder: str, num_workers: int) -> Dict[str, List[str]]:

    # Get all xml paths
    fpath_xmls = []
    for root, dirs, files in os.walk(xml_folder):
        for file in files:
            if os.path.splitext(file)[-1].lower() != XML_EXT:
                continue
            fpath_xmls.append(os.path.join(root, file))
    logging.info(f'Found {len(fpath_xmls)} XMLs at {xml_folder}')

    # Read through xmls to get MRN in parallel
    with multiprocessing.Pool(processes=num_workers) as pool:
        mrn_xml_list = pool.starmap(
            _map_mrn_to_xml,
            [(fpath_xml,) for fpath_xml in fpath_xmls],
        )

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


def _data_from_xml(fpath_xml: str) -> Dict[str, Union[str, Dict[str, np.ndarray]]]:
    ecg_data = dict()

    # define tags that we want to find and use SoupStrainer to speed up search
    tags = [
        'patientdemographics',
        'testdemographics',
        'order',
        'restingecgmeasurements',
        'originalrestingecgmeasurements',
        'diagnosis',
        'originaldiagnosis',
        'intervalmeasurementtimeresolution',
        'intervalmeasurementamplituderesolution',
        'intervalmeasurementfilter',
        'amplitudemeasurements',
        'measurementmatrix',
        'waveform',
    ]
    strainer = bs4.SoupStrainer(tags)

    # lxml parser makes all tags lower case
    with open(fpath_xml, 'r') as f:
        soup = bs4.BeautifulSoup(f, 'lxml', parse_only=strainer)

    for tag in tags:
        tag_suffix = ''
        if tag == 'restingecgmeasurements':
            tag_suffix = '_md'
        elif tag == 'originalrestingecgmeasurements':
            tag_suffix = '_pc'
        elif tag == 'diagnosis':
            soup_tag = soup.find(tag)
            if soup_tag is not None:
                ecg_data['diagnosis_md'] = _parse_soup_diagnosis(soup_tag)
            continue
        elif tag == 'originaldiagnosis':
            soup_tag = soup.find(tag)
            if soup_tag is not None:
                ecg_data['diagnosis_pc'] = _parse_soup_diagnosis(soup_tag)
            continue
        elif tag == 'amplitudemeasurements':
            soup_tag = soup.find(tag)
            if soup_tag is not None:
                amplitude_data = _get_amplitude_from_amplitude_tags(soup.find_all('measuredamplitude'))
                ecg_data['amplitude'] = amplitude_data
            continue
        elif tag == 'measurementmatrix':
            soup_tag = soup.find(tag)
            if soup_tag is not None:
                ecg_data['measurementmatrix'] = _get_measurement_matrix_from_matrix_tags(soup.find_all('measurementmatrix'))
            continue
        elif tag == 'waveform':
            voltage_data = _get_voltage_from_waveform_tags(soup.find_all(tag))
            ecg_data.update(voltage_data)
            continue

        soup_tag = soup.find(tag)

        if soup_tag is not None:
            # find sub tags
            soup_sub_tags = soup_tag.find_all()

            # if there are no sub tags, use original tag
            if len(soup_sub_tags) == 0:
                soup_sub_tags = [soup_tag]

            ecg_data.update({st.name + tag_suffix: st.text for st in soup_sub_tags})

    return ecg_data


def _parse_soup_diagnosis(input_from_soup: bs4.Tag) -> str:

    parsed_text = ''

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


def _get_amplitude_from_amplitude_tags(amplitude_tags: bs4.ResultSet) -> Dict[str, Union[str, Dict[str, np.ndarray]]]:
    wave_ids = set()
    amplitude_data = {}
    amplitude_features = ['peak', 'start', 'duration', 'area']
    ecg_rest_amp_leads = {k.upper(): v for k, v in ECG_REST_AMP_LEADS.items()}
    for amplitude_tag in amplitude_tags:
        lead_id = amplitude_tag.find('amplitudemeasurementleadid').text
        wave_id = amplitude_tag.find('amplitudemeasurementwaveid').text
        if wave_id not in wave_ids:
            wave_ids.add(wave_id)
            for amplitude_feature in amplitude_features:
                amplitude_data[f'measuredamplitude{amplitude_feature}_{wave_id}'] = np.empty(len(ecg_rest_amp_leads))
                amplitude_data[f'measuredamplitude{amplitude_feature}_{wave_id}'][:] = np.nan
        for amplitude_feature in amplitude_features:
            value = int(amplitude_tag.find(f'amplitudemeasurement{amplitude_feature}').text)
            try:
                amplitude_data[f'measuredamplitude{amplitude_feature}_{wave_id}'][ecg_rest_amp_leads[lead_id]] = value
            except KeyError as e:
                logging.warning(f'Amplitude of lead {str(e)} will not be extracted.')
    return amplitude_data


def _decode_array(array_raw: str, scale: float = 1.0) -> np.ndarray:
    decoded = base64.b64decode(array_raw)
    waveform = [struct.unpack("h", bytes([decoded[t], decoded[t + 1]]))[0] for t in range(0, len(decoded), 2)]
    return np.array(waveform) * scale


def _get_measurement_matrix_from_matrix_tags(matrix_tags: bs4.ResultSet) -> np.ndarray:
    for matrix_tag in matrix_tags:
        matrix = matrix_tag.text
        decoded = _decode_array(matrix)
    return decoded


def _get_voltage_from_waveform_tags(waveform_tags: bs4.ResultSet) -> Dict[str, Union[str, Dict[str, np.ndarray]]]:
    voltage_data = dict()
    metadata_tags = ['samplebase', 'sampleexponent', 'highpassfilter', 'lowpassfilter', 'acfilter']

    for waveform_tag in waveform_tags:
        # only use full rhythm waveforms, do not use median waveforms
        if waveform_tag.find('waveformtype').text != 'Rhythm':
            continue

        # get voltage metadata
        for metadata_tag in metadata_tags:
            mt = waveform_tag.find(metadata_tag)
            if mt is not None:
                voltage_data[f'waveform_{metadata_tag}'] = mt.text

        # get voltage leads and lead metadata
        lead_data = _get_voltage_from_lead_tags(waveform_tag.find_all('leaddata'))
        voltage_data.update(lead_data)
        break
    return voltage_data


def _get_voltage_from_lead_tags(lead_tags: bs4.ResultSet) -> Dict[str, Union[str, Dict[str, np.ndarray]]]:
    lead_data = dict()
    voltage = dict()
    all_lead_lengths = []
    all_lead_units = []

    try:
        for lead_tag in lead_tags:
            # for each lead, we make sure all leads use 2 bytes per sample,
            # the decoded lead length is the same as the lead length tag,
            # the lead lengths are all the same, and the units are all the same
            lead_sample_size = int(lead_tag.find('leadsamplesize').text)
            assert lead_sample_size == 2

            lead_id = lead_tag.find('leadid').text
            lead_scale = lead_tag.find('leadamplitudeunitsperbit').text
            lead_waveform_raw = lead_tag.find('waveformdata').text
            lead_waveform = _decode_array(lead_waveform_raw, float(lead_scale))

            lead_length = lead_tag.find('leadsamplecounttotal').text
            lead_units = lead_tag.find('leadamplitudeunits').text

            assert int(lead_length) == len(lead_waveform)
            all_lead_lengths.append(lead_length)
            all_lead_units.append(lead_units)

            voltage[lead_id] = lead_waveform

        # vector math to get remaining leads
        assert len(voltage) == 8
        voltage['III'] = voltage['II'] - voltage['I']
        voltage['aVR'] = -1 * (voltage['I'] + voltage['II']) / 2
        voltage['aVL'] = voltage['I'] - voltage['II'] / 2
        voltage['aVF'] = voltage['II'] - voltage['I'] / 2

        # add voltage length and units to metadata
        assert len(set(all_lead_lengths)) == 1
        assert len(set(all_lead_units)) == 1
        lead_data['voltagelength'] = all_lead_lengths[0]
        lead_data['voltageunits'] = all_lead_units[0]

        lead_data['voltage'] = voltage
        return lead_data
    except (AssertionError, AttributeError, ValueError) as e:
        logging.exception(e)
        return dict()


def _compress_and_save_data(
    hd5: h5py.Group, name: str, data: Union[str, np.ndarray],
    dtype: str, method: str = 'zstd', compression_opts: int = 19,
) -> None:
    # Define codec
    codec = numcodecs.zstd.Zstd(level=compression_opts)

    # If data is string, encode to bytes
    if dtype == 'str':
        # do not save empty string, cannot be decoded
        if data == '':
            return
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

    # Extract data from XML into dict
    ecg_data = _data_from_xml(fpath_xml)
    dt = datetime.strptime(f"{ecg_data['acquisitiondate']} {ecg_data['acquisitiontime']}", '%m-%d-%Y %H:%M:%S')
    ecg_dt = dt.isoformat()

    if (os.stat(fpath_xml).st_size == 0 or not ecg_data):
        # If XML is empty, remove the XML file and do not convert
        os.remove(fpath_xml)
        convert = 0
        logging.warning(f'Conversion of {fpath_xml} failed! XML is empty.')
    elif ecg_dt in hd5.keys():
        # If patient already has an ECG at given date and time, skip duplicate
        logging.warning(f'Conversion of {fpath_xml} skipped. Converted XML already exists in HD5.')
        convert = -1
    elif "voltage" not in ecg_data:
        # If we could not get voltage, do not convert (see _get_voltage_from_lead_tags)
        logging.warning(f'Conversion of {fpath_xml} failed! Voltage is empty or badly formatted.')
        convert = 0
    elif _get_max_voltage(ecg_data["voltage"]) == 0:
        # If the max voltage value is 0, do not convert
        logging.warning(f'Conversion of {fpath_xml} failed! Maximum voltage is 0.')
        convert = 0

    # If all prior checks passed, write hd5 group for ECG
    if convert == 1:
        gp = hd5.create_group(ecg_dt)

        # Save voltage leads
        voltage = ecg_data.pop('voltage')
        for lead in voltage:
            _compress_and_save_data(hd5=gp, name=lead, data=voltage[lead].astype('int16'), dtype='int16')

        # Save ECG wave amplitudes if present
        try:
            amplitudes = ecg_data.pop('amplitude')
            for amplitude in amplitudes:
                _compress_and_save_data(hd5=gp, name=amplitude, data=amplitudes[amplitude], dtype='float')
        except KeyError:
            logging.info(f'Conversion of amplitude measures failed! Not present in: {fpath_xml}')
        
        # Save measurement matrix if present
        try:
            measurement_matrix = ecg_data.pop('measurementmatrix')
            _compress_and_save_data(hd5=gp, name='measurementmatrix', data=measurement_matrix.astype('int16'), dtype='int16')
        except KeyError:
            logging.info(f'Conversion of measurement matrix failed! Not present in: {fpath_xml}')

        # Save everything else
        for key in ecg_data:
            _compress_and_save_data(hd5=gp, name=key, data=ecg_data[key], dtype='str')

        # Clean Patient MRN to only numbers
        key_mrn_clean = 'patientid_clean'
        if 'patientid' in ecg_data:
            mrn_clean = _clean_mrn(ecg_data['patientid'], fallback='')
            _compress_and_save_data(hd5=gp, name=key_mrn_clean, data=mrn_clean, dtype='str')

        # Clean cardiologist read
        key_read_md = 'diagnosis_md'
        key_read_md_clean = 'read_md_clean'
        if key_read_md in ecg_data:
            read_md_clean = _clean_read_text(text=ecg_data[key_read_md])
            _compress_and_save_data(hd5=gp, name=key_read_md_clean, data=read_md_clean, dtype='str')

        # Clean MUSE read
        key_read_pc = 'diagnosis_pc'
        key_read_pc_clean = 'read_pc_clean'
        if key_read_pc in ecg_data:
            read_pc_clean = _clean_read_text(text=ecg_data[key_read_pc])
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


def _convert_mrn_xmls_to_hd5_wrapper(mrn_xmls_map: Dict[str, List[str]], dir_hd5: str, num_workers: int, hd5_prefix: str = 'partners_ecg_rest'):
    tot_xml = sum([len(v) for k, v in mrn_xmls_map.items()])
    os.makedirs(dir_hd5, exist_ok=True)

    with multiprocessing.Pool(processes=num_workers) as pool:
        converted = pool.starmap(
            _convert_mrn_xmls_to_hd5,
            [(mrn, fpath_xmls, dir_hd5, hd5_prefix) for mrn, fpath_xmls in mrn_xmls_map.items()],
        )
    num_hd5 = sum([x[0] for x in converted])
    num_xml = sum([x[1] for x in converted])
    num_dup = sum([x[2] for x in converted])

    logging.info(f"Converted {num_xml} XMLs to {num_hd5} HD5s at {dir_hd5}")
    logging.info(f"Skipped {num_dup} duplicate XMLs")
    logging.info(f"Skipped {tot_xml - num_dup - num_xml} malformed XMLs")
