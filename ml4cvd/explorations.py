# explorations.py

# Imports
import os
import csv
import math
import operator
import datetime
from functools import reduce
from itertools import combinations
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Generator, Optional, DefaultDict

import h5py
import logging
import numpy as np
import pandas as pd
from tensorflow.keras.models import Model

import matplotlib
matplotlib.use('Agg')  # Need this to write images from the GSA servers.  Order matters:
import matplotlib.pyplot as plt  # First import matplotlib, then use Agg, then import plt

from ml4cvd.TensorMap import TensorMap
from ml4cvd.tensor_generators import TensorGenerator
from ml4cvd.models import make_multimodal_multitask_model
from ml4cvd.plots import plot_histograms_in_pdf, plot_heatmap, evaluate_predictions, subplot_rocs, subplot_scatters
from ml4cvd.defines import TENSOR_EXT, IMAGE_EXT, ECG_CHAR_2_IDX, ECG_IDX_2_CHAR, CODING_VALUES_MISSING, CODING_VALUES_LESS_THAN_ONE, JOIN_CHAR, MRI_SEGMENTED_CHANNEL_MAP

CSV_EXT = '.tsv'


def sort_csv(input_csv_file, value_csv):
    lvef = {}
    with open(value_csv, 'r') as value_file:
        lol = list(csv.reader(value_file, delimiter='\t'))
        logging.info('CSV of MRI volumes header:{}'.format(list(enumerate(lol[0]))))
        for row in lol[1:]:
            sample_id = row[0]
            if row[5] != 'NA':
                lvef[sample_id] = float(row[5])

    print('try:', input_csv_file.replace(CSV_EXT, '_diff_sorted'+CSV_EXT))
    with open(input_csv_file, mode='r') as input_csv:
        with open(input_csv_file.replace(CSV_EXT, '_diff_sorted'+CSV_EXT), mode='w') as output_csv:
            csv_writer = csv.writer(output_csv, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            csv_reader = csv.reader(input_csv, delimiter='\t')
            csv_writer.writerow(next(csv_reader)+['discrepancy'])
            csv_sorted = sorted(csv_reader, key=lambda row: abs(float(lvef[row[0]])-float(row[5])), reverse=True)
            [csv_writer.writerow(row + [float(lvef[row[0]])-float(row[5])]) for row in csv_sorted]


def predictions_to_pngs(predictions: np.ndarray, tensor_maps_in: List[TensorMap], tensor_maps_out: List[TensorMap], data: Dict[str, np.ndarray],
                        labels: Dict[str, np.ndarray], paths: List[str], folder: str) -> None:
    # TODO Remove this command line order dependency
    input_map = tensor_maps_in[0]
    if not os.path.exists(folder):
        os.makedirs(folder)
    for y, tm in zip(predictions, tensor_maps_out):
        if not isinstance(predictions, list):  # When models have a single output model.predict returns a ndarray otherwise it returns a list
            y = predictions
        for im in tensor_maps_in:
            if tm.is_categorical() and im.dependent_map == tm:
                input_map = im
            elif len(tm.shape) == len(im.shape):
                input_map = im
        logging.info(f"Write predictions as PNGs y:{y.shape} labels:{labels[tm.output_name()].shape} folder:{folder}")
        if tm.is_mesh():
            vmin = np.min(data[input_map.input_name()])
            vmax = np.max(data[input_map.input_name()])
            for i in range(y.shape[0]):
                sample_id = os.path.basename(paths[i]).replace(TENSOR_EXT, '')
                if input_map.axes() == 4 and input_map.shape[-1] == 1:
                    sample_data = data[input_map.input_name()][i, ..., 0]
                    cols = max(2, int(math.ceil(math.sqrt(sample_data.shape[-1]))))
                    rows = max(2, int(math.ceil(sample_data.shape[-1] / cols)))
                    path_prefix = f'{folder}{sample_id}_bbox_batch_{i:02d}{IMAGE_EXT}'
                    logging.info(f"sample_data shape: {sample_data.shape} cols {cols}, {rows} Predicted BBox: {y[i]}, True BBox: {labels[tm.output_name()][i]} Vmin {vmin} Vmax{vmax}")
                    _plot_3d_tensor_slices_as_gray(sample_data, path_prefix, cols, rows, bboxes=[labels[tm.output_name()][i], y[i]])
                else:
                    fig, ax = plt.subplots(1)
                    if input_map.axes() == 3 and input_map.shape[-1] == 1:
                        ax.imshow(data[input_map.input_name()][i, :, :, 0], cmap='gray', vmin=vmin, vmax=vmax)
                    elif input_map.axes() == 2:
                        ax.imshow(data[input_map.input_name()][i, :, :], cmap='gray', vmin=vmin, vmax=vmax)
                    corner, width, height = _2d_bbox_to_corner_and_size(labels[tm.output_name()][i])
                    ax.add_patch(matplotlib.patches.Rectangle(corner, width, height, linewidth=1, edgecolor='g', facecolor='none'))
                    y_corner, y_width, y_height = _2d_bbox_to_corner_and_size(y[i])
                    ax.add_patch(matplotlib.patches.Rectangle(y_corner, y_width, y_height, linewidth=1, edgecolor='y', facecolor='none'))
                    logging.info(f"True BBox: {corner}, {width}, {height} Predicted BBox: {y_corner}, {y_width}, {y_height} Vmin {vmin} Vmax{vmax}")
                plt.savefig(f"{folder}{sample_id}_bbox_batch_{i:02d}{IMAGE_EXT}")
        elif len(tm.shape) in [1, 2]:
            vmin = np.min(data[input_map.input_name()])
            vmax = np.max(data[input_map.input_name()])
            for i in range(y.shape[0]):
                sample_id = os.path.basename(paths[i]).replace(TENSOR_EXT, '')
                if input_map.axes() == 3 and input_map.shape[-1] == 1:
                    plt.imsave(f"{folder}{sample_id}_batch_{i:02d}{IMAGE_EXT}", data[input_map.input_name()][i, :, :, 0], cmap='gray', vmin=vmin, vmax=vmax)
                elif input_map.axes() == 2:
                    plt.imsave(f"{folder}{sample_id}_batch_{i:02d}{IMAGE_EXT}", data[input_map.input_name()][i, :, :], cmap='gray', vmin=vmin, vmax=vmax)
                elif input_map.axes() == 3:
                    for j in range(data[input_map.input_name()].shape[-1]):
                        image_file = f"{folder}{sample_id}_batch_{i:02d}_slice_{j:02d}{IMAGE_EXT}"
                        plt.imsave(image_file, data[input_map.input_name()][i, :, :, j], cmap='gray', vmin=vmin, vmax=vmax)
                elif input_map.axes() == 4:
                    for j in range(data[input_map.input_name()].shape[-2]):
                        image_file = f"{folder}{sample_id}_batch_{i:02d}_slice_{j:02d}{IMAGE_EXT}"
                        plt.imsave(image_file, data[input_map.input_name()][i, :, :, j, 0], cmap='gray', vmin=vmin, vmax=vmax)
        elif len(tm.shape) == 3:
            for i in range(y.shape[0]):
                sample_id = os.path.basename(paths[i]).replace(TENSOR_EXT, '')
                if tm.is_categorical():
                    plt.imsave(f"{folder}{sample_id}_truth_{i:02d}{IMAGE_EXT}", np.argmax(labels[tm.output_name()][i], axis=-1, cmap='gray'))
                    plt.imsave(f"{folder}{sample_id}_prediction_{i:02d}{IMAGE_EXT}", np.argmax(y[i], axis=-1), cmap='gray')
                    if input_map is not None:
                        plt.imsave(f"{folder}{sample_id}_mri_slice_{i:02d}{IMAGE_EXT}", data[input_map.input_name()][i, :, :, 0], cmap='gray')
                else:
                    for j in range(y.shape[3]):
                        plt.imsave(f"{folder}{sample_id}_truth_{i:02d}_{j:02d}{IMAGE_EXT}", labels[tm.output_name()][i, :, :, j], cmap='gray')
                        plt.imsave(f"{folder}{sample_id}_prediction_{i:02d}_{j:02d}{IMAGE_EXT}", y[i, :, :, j], cmap='gray')
                        plt.imsave(f"{folder}{sample_id}_mri_slice_{i:02d}_{j:02d}{IMAGE_EXT}", data[input_map.input_name()][i, :, :, j], cmap='gray')
        elif len(tm.shape) == 4:
            for im in tensor_maps_in:
                if im.dependent_map == tm:
                    break
            for i in range(y.shape[0]):
                sample_id = os.path.basename(paths[i]).replace(TENSOR_EXT, '')
                for j in range(y.shape[3]):
                    if tm.is_categorical():
                        truth = np.argmax(labels[tm.output_name()][i, :, :, j, :], axis=-1)
                        prediction = np.argmax(y[i, :, :, j, :], axis=-1)
                        true_donut = np.ma.masked_where(truth == 2, data[im.input_name()][i, :, :, j, 0])
                        predict_donut = np.ma.masked_where(prediction == 2, data[im.input_name()][i, :, :, j, 0])
                        plt.imsave(folder+sample_id+'_truth_{0:03d}_{1:03d}'.format(i, j)+IMAGE_EXT, truth, cmap='gray')
                        plt.imsave(folder+sample_id+'_prediction_{0:03d}_{1:03d}'.format(i, j)+IMAGE_EXT, prediction, cmap='gray')
                        plt.imsave(folder+sample_id+'_mri_slice_{0:03d}_{1:03d}'.format(i, j)+IMAGE_EXT, data[im.input_name()][i, :, :, j, 0], cmap='gray')
                        plt.imsave(folder+sample_id + '_true_donut_{0:03d}_{1:03d}'.format(i, j) + IMAGE_EXT, true_donut, cmap='gray')
                        plt.imsave(folder + sample_id + '_predict_donut_{0:03d}_{1:03d}'.format(i, j) + IMAGE_EXT, predict_donut, cmap='gray')
                    else:
                        plt.imsave(folder+sample_id+'_truth_{0:03d}_{1:03d}'.format(i, j)+IMAGE_EXT, labels[tm.output_name()][i, :, :, j, 0], cmap='gray')
                        plt.imsave(folder+sample_id+'_prediction_{0:03d}_{1:03d}'.format(i, j)+IMAGE_EXT, y[i, :, :, j, 0], cmap='gray')


def plot_while_learning(model, tensor_maps_in: List[TensorMap], tensor_maps_out: List[TensorMap],
                        generate_train: Generator[Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Optional[List[str]]], None, None],
                        test_data: Dict[str, np.ndarray], test_labels: Dict[str, np.ndarray], test_paths: List[str], epochs: int, batch_size: int,
                        training_steps: int, folder: str, write_pngs: bool):
    if not os.path.exists(folder):
        os.makedirs(folder)

    for i in range(epochs):
        rocs = []
        scatters = []
        predictions = model.predict(test_data, batch_size=batch_size)
        if len(tensor_maps_out) == 1:
            predictions = [predictions]
        for y, tm in zip(predictions, tensor_maps_out):
            for im in tensor_maps_in:
                if im.dependent_map == tm:
                    break
            if not write_pngs:
                mri_in = test_data[im.input_name()]
                vmin = np.min(mri_in)
                vmax = np.max(mri_in)
                logging.info(f"epoch:{i} write segmented mris y shape:{y.shape} label shape:{test_labels[tm.output_name()].shape} to folder:{folder}")
                if tm.is_categorical() and len(tm.shape) == 3:
                    for yi in range(y.shape[0]):
                        plt.imsave(f"{folder}batch_{yi}_truth_epoch_{i:03d}{IMAGE_EXT}", np.argmax(test_labels[tm.output_name()][yi], axis=-1), cmap='gray')
                        plt.imsave(f"{folder}batch_{yi}_prediction_epoch_{i:03d}{IMAGE_EXT}", np.argmax(y[yi], axis=-1), cmap='gray')
                        plt.imsave(f"{folder}batch_{yi}_mri_epoch_{i:03d}{IMAGE_EXT}", mri_in[yi, :, :, 0], cmap='gray', vmin=vmin, vmax=vmax)
                elif tm.is_categorical() and len(tm.shape) == 4:
                    for yi in range(y.shape[0]):
                        for j in range(y.shape[3]):
                            truth = np.argmax(test_labels[tm.output_name()][yi, :, :, j, :], axis=-1)
                            prediction = np.argmax(y[yi, :, :, j, :], axis=-1)
                            true_donut = np.ma.masked_where(truth == 2, mri_in[yi, :, :, j, 0])
                            predict_donut = np.ma.masked_where(prediction == 2, mri_in[yi, :, :, j, 0])
                            plt.imsave(f"{folder}batch_{yi}_slice_{j:03d}_prediction_epoch_{i:03d}{IMAGE_EXT}", prediction, cmap='gray')
                            plt.imsave(f"{folder}batch_{yi}_slice_{j:03d}_p_donut_epoch_{i:03d}{IMAGE_EXT}", predict_donut, cmap='gray', vmin=vmin, vmax=vmax)
                            if i == 0:
                                plt.imsave(f"{folder}batch_{yi}_slice_{j:03d}_truth_epoch_{i:03d}{IMAGE_EXT}", truth, cmap='gray')
                                plt.imsave(f"{folder}batch_{yi}_slice_{j:03d}_t_donut_epoch_{i:03d}{IMAGE_EXT}", true_donut, cmap='gray', vmin=vmin, vmax=vmax)
                                plt.imsave(f"{folder}batch_{yi}_slice_{j:03d}_mri_epoch_{i:03d}{IMAGE_EXT}", mri_in[yi, :, :, j, 0], cmap='gray', vmin=vmin, vmax=vmax)
                else:
                    logging.warning(f'Not writing PNGs')
            elif write_pngs:
                if len(tensor_maps_out) == 1:
                    y = predictions[0]
                evaluate_predictions(tm, y, test_labels[tm.output_name()], f"{tm.name}_epoch_{i:03d}", folder, test_paths, rocs=rocs, scatters=scatters)
        if len(rocs) > 1:
            subplot_rocs(rocs, folder+f"epoch_{i:03d}_")
        if len(scatters) > 1:
            subplot_scatters(scatters, folder+f"epoch_{i:03d}_")

        model.fit_generator(generate_train, steps_per_epoch=training_steps, epochs=1, verbose=1)


def plot_histograms_of_tensors_in_pdf(run_id: str,
                                      tensor_folder: str,
                                      output_folder: str,
                                      max_samples: int = None) -> None:
    """
    :param id: name for the plotting run
    :param tensor_folder: directory with tensor files to plot histograms from
    :param output_folder: folder containing the output plot
    :param max_samples: specifies how many tensor files to down-sample from; by default all tensors are used
    """
    stats, num_tensor_files = _collect_continuous_stats_from_tensor_files(tensor_folder, max_samples)
    logging.info(f"Collected continuous stats for {len(stats)} fields. Now plotting histograms of them...")
    plot_histograms_in_pdf(stats, num_tensor_files, run_id, output_folder)


def plot_heatmap_of_tensors(id: str,
                            tensor_folder: str,
                            output_folder: str,
                            min_samples: int,
                            max_samples: int = None) -> None:
    """
    :param id: name for the plotting run
    :param tensor_folder: directory with tensor files to plot histograms from
    :param output_folder: folder containing the output plot
    :param min_samples: calculate correlation coefficient only if both fields have values from that many common samples
    :param max_samples: specifies how many tensor files to down-sample from; by default all tensors are used
    """
    stats, _ = _collect_continuous_stats_from_tensor_files(tensor_folder, max_samples, ['0'], 0)
    logging.info(f"Collected continuous stats for {len(stats)} fields. Now plotting a heatmap of their correlations...")
    plot_heatmap(stats, id, min_samples, output_folder)


def tabulate_correlations_of_tensors(run_id: str,
                                     tensor_folder: str,
                                     output_folder: str,
                                     min_samples: int,
                                     max_samples: int = None) -> None:
    """
    :param id: name for the plotting run
    :param tensor_folder: directory with tensor files to plot histograms from
    :param output_folder: folder containing the output plot
    :param min_samples: calculate correlation coefficient only if both fields have values from that many common samples
    :param max_samples: specifies how many tensor files to down-sample from; by default all tensors are used
    """
    stats, _ = _collect_continuous_stats_from_tensor_files(tensor_folder, max_samples)
    logging.info(f"Collected continuous stats for {len(stats)} fields. Now tabulating their cross-correlations...")
    _tabulate_correlations(stats, run_id, min_samples, output_folder)


def mri_dates(tensors: str, output_folder: str, run_id: str):
    incident_dates = []
    prevalent_dates = []
    disease = 'hypertension'
    disease_date_key = disease + '_date'
    data_date_key = 'assessment-date_0_0'
    tensor_paths = [tensors + tp for tp in os.listdir(tensors) if os.path.splitext(tp)[-1].lower() == TENSOR_EXT]
    for tp in tensor_paths:
        try:
            with h5py.File(tp, 'r') as hd5:
                if data_date_key in hd5 and disease_date_key in hd5:
                    if int(hd5[disease][0]) == 1:
                        data_date = str2date(str(hd5[data_date_key][0]))
                        disease_date = str2date(str(hd5[disease_date_key][0]))
                        if data_date < disease_date:
                            incident_dates.append(disease_date)
                        else:
                            prevalent_dates.append(disease_date)
        except:
            logging.exception(f"Broken tensor at:{tp}")

    plt.figure(figsize=(12, 12))
    plt.xlabel(data_date_key)
    plt.hist(incident_dates, bins=60)
    plt.savefig(os.path.join(output_folder, run_id, disease+'_'+data_date_key+'_incident'+IMAGE_EXT))
    plt.figure(figsize=(12,12))
    plt.xlabel(data_date_key)
    plt.hist(prevalent_dates, bins=60)
    plt.savefig(os.path.join(output_folder, run_id, disease+'_'+data_date_key+'_prevalent'+IMAGE_EXT))


def ecg_dates(tensors: str, output_folder: str, run_id: str):
    incident_dates = []
    prevalent_dates = []
    tensor_paths = [tensors + tp for tp in os.listdir(tensors) if os.path.splitext(tp)[-1].lower()==TENSOR_EXT]
    for tp in tensor_paths:
        try:
            with h5py.File(tp, 'r') as hd5:
                if 'ecg_bike_date_0' in hd5 and 'coronary_artery_disease_soft_date' in hd5:
                    ecg_date = str2date(str(hd5['ecg_bike_date_0'][0]))
                    cad_date = str2date(str(hd5['coronary_artery_disease_soft_date'][0]))
                    if ecg_date < cad_date:
                        incident_dates.append(ecg_date)
                    else:
                        prevalent_dates.append(ecg_date)
        except:
            logging.exception(f"Broken tensor at:{tp}")

    plt.figure(figsize=(12, 12))
    plt.xlabel('ECG Acquisition Date')
    plt.hist(incident_dates, bins=60)
    plt.savefig(os.path.join(output_folder, run_id, 'ecg_dates_incident'+IMAGE_EXT))
    plt.figure(figsize=(12, 12))
    plt.xlabel('ECG Acquisition Date')
    plt.hist(prevalent_dates, bins=60)
    plt.savefig(os.path.join(output_folder, run_id, 'ecg_dates_prevalent'+IMAGE_EXT))


def str2date(d):
    parts = d.split('-')
    if len(parts) < 2:
        return datetime.datetime.now().date()
    return datetime.date(int(parts[0]), int(parts[1]), int(parts[2]))


def sample_from_char_model(char_model: Model, test_batch: Dict[str, np.ndarray], test_paths: List[str]) -> None:
    window_size = test_batch['input_ecg_rest_text_ecg_text'].shape[1]
    alphabet_size = test_batch['input_ecg_rest_text_ecg_text'].shape[2]
    for i in range(test_batch['input_embed_hidden_layer'].shape[0]):
        count = 0
        sentence = ''
        next_char = ''
        embed_in = test_batch['input_embed_hidden_layer'][i:i+1, :]
        burn_in = np.zeros((1, window_size, alphabet_size), dtype=np.float32)
        window_size = burn_in.shape[1]
        with h5py.File(test_paths[i], 'r') as hd5:
            logging.info(f"\n")
            logging.info(f"Real text: {str(hd5['ecg_rest_text'][0]).strip()}")
        while next_char != '!' and count < 400:
            cur_test = {'input_embed_hidden_layer': embed_in, 'input_ecg_rest_text_ecg_text': burn_in}
            y_pred = char_model.predict(cur_test)
            next_char = ECG_IDX_2_CHAR[_sample_with_heat(y_pred[0, :], 0.7)]
            sentence += next_char
            burn_in = np.zeros((1,) + test_batch['input_ecg_rest_text_ecg_text'].shape[1:], dtype=np.float32)
            for j, c in enumerate(reversed(sentence)):
                if j == window_size:
                    break
                burn_in[0, window_size-j-1, ECG_CHAR_2_IDX[c]] = 1.0
            count += 1
        logging.info(f"Model text:{sentence}")


def tensors_to_label_dictionary(categorical_labels: List,
                                continuous_labels: List,
                                gene_labels: List,
                                samples2genes: Dict[str, str],
                                test_paths: List) -> Dict[str, np.ndarray]:
    label_dict = {k: np.zeros((len(test_paths))) for k in categorical_labels + continuous_labels + gene_labels}
    for i, tp in enumerate(test_paths):
        hd5 = h5py.File(tp, 'r')
        for k in categorical_labels:
            if k in hd5['categorical']:
                label_dict[k][i] = 1
            elif k in hd5 and hd5[k][0] == 1:
                label_dict[k][i] = 1
        for mk in continuous_labels:
            for k in mk.split('|'):
                if k in hd5['continuous']:
                    label_dict[mk][i] = hd5['continuous'][k][0]
        for k in gene_labels:
            if tp in samples2genes and samples2genes[tp] == k:
                label_dict[k][i] = 1

    return label_dict


def test_labels_to_label_map(test_labels: Dict[TensorMap, np.ndarray], examples: int) -> Tuple[Dict[str, np.ndarray], List[str], List[str]]:
    label_dict = {tm: np.zeros((examples,)) for tm in test_labels}
    categorical_labels = []
    continuous_labels = []

    for tm in test_labels:
        for i in range(examples):
            if tm.is_continuous() and tm.axes() == 1:
                label_dict[tm][i] = tm.rescale(test_labels[tm][i])
                continuous_labels.append(tm)
            elif tm.is_categorical() and tm.axes() == 1:
                label_dict[tm][i] = np.argmax(test_labels[tm][i])
                categorical_labels.append(tm)

    return label_dict, categorical_labels, continuous_labels


def infer_with_pixels(args):
    stats = Counter()
    tensor_paths_inferred = {}
    args.num_workers = 0
    inference_tsv = os.path.join(args.output_folder, args.id, 'pixel_inference_' + args.id + '.tsv')
    tensor_paths = [args.tensors + tp for tp in sorted(os.listdir(args.tensors)) if os.path.splitext(tp)[-1].lower() == TENSOR_EXT]
    # hard code batch size to 1 so we can iterate over file names and generated tensors together in the tensor_paths for loop
    model = make_multimodal_multitask_model(**args.__dict__)
    generate_test = TensorGenerator(1, args.tensor_maps_in, args.tensor_maps_out, tensor_paths, num_workers=args.num_workers,
                                    cache_size=args.cache_size, keep_paths=True, mixup=args.mixup_alpha)
    with open(inference_tsv, mode='w') as inference_file:
        inference_writer = csv.writer(inference_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        header = ['sample_id']
        for ot, otm in zip(args.output_tensors, args.tensor_maps_out):
            if len(otm.shape) == 1 and otm.is_continuous():
                header.extend([ot+'_prediction', ot+'_actual'])
            elif len(otm.shape) == 1 and otm.is_categorical():
                channel_columns = []
                for k in otm.channel_map:
                    channel_columns.append(ot + '_' + k + '_prediction')
                    channel_columns.append(ot + '_' + k + '_actual')
                header.extend(channel_columns)
            elif otm.name in ['mri_systole_diastole_8_segmented', 'sax_all_diastole_segmented']:
                pix_tm = args.tensor_maps_in[1]
                header.extend(['pixel_size', 'background_pixel_prediction', 'background_pixel_actual', 'ventricle_pixel_prediction', 'ventricle_pixel_actual', 'myocardium_pixel_prediction', 'myocardium_pixel_actual'])
                if otm.name == 'sax_all_diastole_segmented':
                    header.append('total_b_slices')
        inference_writer.writerow(header)

        while True:
            input_data, true_label, tensor_path = next(generate_test)
            if tensor_path[0] in tensor_paths_inferred:
                logging.info(f"Inference on {stats['count']} tensors finished. Inference TSV file at: {inference_tsv}")
                break

            prediction = model.predict(input_data)
            if len(args.tensor_maps_out) == 1:
                prediction = [prediction]

            csv_row = [os.path.basename(tensor_path[0]).replace(TENSOR_EXT, '')]  # extract sample id
            for y, tm in zip(prediction, args.tensor_maps_out):
                if len(tm.shape) == 1 and tm.is_continuous():
                    csv_row.append(str(tm.rescale(y)[0][0]))  # first index into batch then index into the 1x1 structure
                    if tm.sentinel is not None and tm.sentinel == true_label[tm.output_name()][0][0]:
                        csv_row.append("NA")
                    else:
                        csv_row.append(str(tm.rescale(true_label[tm.output_name()])[0][0]))
                elif len(tm.shape) == 1 and tm.is_categorical():
                    for k in tm.channel_map:
                        csv_row.append(str(y[0][tm.channel_map[k]]))
                        csv_row.append(str(true_label[tm.output_name()][0][tm.channel_map[k]]))
                elif tm.name in ['mri_systole_diastole_8_segmented', 'sax_all_diastole_segmented']:
                    csv_row.append(f"{pix_tm.rescale(input_data['input_mri_pixel_width_cine_segmented_sax_inlinevf_continuous'][0][0]):0.3f}")
                    csv_row.append(f'{np.sum(np.argmax(y, axis=-1) == MRI_SEGMENTED_CHANNEL_MAP["background"]):0.2f}')
                    csv_row.append(f'{np.sum(true_label[tm.output_name()][..., MRI_SEGMENTED_CHANNEL_MAP["background"]]):0.1f}')
                    csv_row.append(f'{np.sum(np.argmax(y, axis=-1) == MRI_SEGMENTED_CHANNEL_MAP["ventricle"]):0.2f}')
                    csv_row.append(f'{np.sum(true_label[tm.output_name()][..., MRI_SEGMENTED_CHANNEL_MAP["ventricle"]]):0.1f}')
                    csv_row.append(f'{np.sum(np.argmax(y, axis=-1) == MRI_SEGMENTED_CHANNEL_MAP["myocardium"]):0.2f}')
                    csv_row.append(f'{np.sum(true_label[tm.output_name()][..., MRI_SEGMENTED_CHANNEL_MAP["myocardium"]]):0.1f}')
                    if tm.name == 'sax_all_diastole_segmented':
                        background_counts = np.count_nonzero(true_label[tm.output_name()][..., MRI_SEGMENTED_CHANNEL_MAP["background"]] == 0, axis=(0, 1, 2))
                        csv_row.append(f'{np.count_nonzero(background_counts):0.0f}')

            inference_writer.writerow(csv_row)
            tensor_paths_inferred[tensor_path[0]] = True
            stats['count'] += 1
            if stats['count'] % 250 == 0:
                logging.info(f"Wrote:{stats['count']} rows of inference.  Last tensor:{tensor_path[0]}")


def _sample_with_heat(preds, temperature=1.0):
    # helper function to sample an index from a probability array
    preds = np.asarray(preds).astype('float64')
    preds = np.log(preds) / temperature
    exp_preds = np.exp(preds)
    preds = exp_preds / np.sum(exp_preds)
    probas = np.random.multinomial(1, preds, 1)
    return np.argmax(probas)


def _2d_bbox_to_corner_and_size(bbox):
    total_axes = bbox.shape[-1] // 2
    lower_left_corner = (bbox[1], bbox[0])
    height = bbox[total_axes] - bbox[0]
    width = bbox[total_axes+1] - bbox[1]
    return lower_left_corner, width, height


def _plot_3d_tensor_slices_as_gray(tensor, figure_path, cols=3, rows=10, bboxes=[]):
    colors = ['blue', 'red', 'green', 'yellow']
    _, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4))
    vmin = np.min(tensor)
    vmax = np.max(tensor)
    for i in range(tensor.shape[-1]):
        axes[i // cols, i % cols].imshow(tensor[:, :, i], cmap='gray', vmin=vmin, vmax=vmax)
        axes[i // cols, i % cols].set_yticklabels([])
        axes[i // cols, i % cols].set_xticklabels([])
        for c, bbox in enumerate(bboxes):
            corner, width, height = _2d_bbox_to_corner_and_size(bbox)
            axes[i // cols, i % cols].add_patch(matplotlib.patches.Rectangle(corner, width, height, linewidth=1, edgecolor=colors[c], facecolor='none'))

    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)


def _tabulate_correlations(stats: Dict[str, Dict[str, List[float]]],
                           output_file_name: str,
                           min_samples: int,
                           output_folder_path: str) -> None:

    """
    Tabulate in pdf correlations of field values given in 'stats'
    :param stats: field names extracted from hd5 dataset names to list of values, one per sample_instance_arrayidx
    :param output_file_name: name of output file in pdf
    :param output_folder_path: directory that output file will be written to
    :param min_samples: calculate correlation coefficient only if both fields have values from that many common samples
    :return: None
    """

    fields = stats.keys()
    num_fields = len(fields)
    field_pairs = combinations(fields, 2)
    table_rows: List[list] = []
    logging.info(f"There are {int(num_fields * (num_fields - 1) / 2)} field pairs.")
    processed_field_pair_count = 0
    nan_counter = Counter()  # keep track of if we've seen a field have NaNs
    for field1, field2 in field_pairs:
        if field1 not in nan_counter.keys() and field2 not in nan_counter.keys():
            common_samples = set(stats[field1].keys()).intersection(stats[field2].keys())
            num_common_samples = len(common_samples)
            processed_field_pair_count += 1
            if processed_field_pair_count % 50000 == 0:
                logging.debug(f"Processed {processed_field_pair_count} field pairs.")
            if num_common_samples >= min_samples:
                field1_values = reduce(operator.concat, [stats[field1][sample] for sample in common_samples])
                field2_values = reduce(operator.concat, [stats[field2][sample] for sample in common_samples])

                num_field1_nans = len(list(filter(math.isnan, field1_values)))
                num_field2_nans = len(list(filter(math.isnan, field2_values)))
                at_least_one_field_has_nans = False
                if num_field1_nans != 0:
                    nan_counter[field1] = True
                    at_least_one_field_has_nans = True
                if num_field2_nans != 0:
                    nan_counter[field2] = True
                    at_least_one_field_has_nans = True
                if at_least_one_field_has_nans:
                    continue

                if len(field1_values) == len(field2_values):
                    if len(set(field1_values)) == 1 or len(set(field2_values)) == 1:
                        logging.debug(f"Not calculating correlation for fields {field1} and {field2} because at least one of "
                                      f"the fields has all the same values for the {num_common_samples} common samples.")
                        continue
                    corr = np.corrcoef(field1_values, field2_values)[1, 0]
                    if not math.isnan(corr):
                        table_rows.append([field1, field2, corr, corr * corr, num_common_samples])
                    else:
                        logging.warning(f"Pearson correlation for fields {field1} and {field2} is NaN.")
                else:
                    logging.debug(f"Not calculating correlation for fields '{field1}' and '{field2}' "
                                  f"because they have different number of values ({len(field1_values)} vs. {len(field2_values)}).")
        else:
            continue

    # Note: NaNs mess up sorting unless they are handled specially by a custom sorting function
    sorted_table_rows = sorted(table_rows, key=operator.itemgetter(2), reverse=True)
    logging.info(f"Total number of correlations: {len(sorted_table_rows)}")

    fields_with_nans = nan_counter.keys()
    if len(fields_with_nans) != 0:
        logging.warning(f"The {len(fields_with_nans)} fields containing NaNs are: {', '.join(fields_with_nans)}.")

    table_path = os.path.join(output_folder_path, output_file_name + CSV_EXT)
    table_header = ["Field 1", "Field 2", "Pearson R", "Pearson R^2",  "Sample Size"]
    df = pd.DataFrame(sorted_table_rows, columns=table_header)
    df.to_csv(table_path, index=False)

    logging.info(f"Saved correlations table at: {table_path}")


def _collect_continuous_stats_from_tensor_files(tensor_folder: str,
                                                max_samples: int = None,
                                                instances: List[str] = ['0', '1', '2'],
                                                max_arr_idx: int = None) -> Tuple[DefaultDict[str, DefaultDict[str, List[float]]], int]:
    if not os.path.exists(tensor_folder):
        raise ValueError('Source directory does not exist: ', tensor_folder)
    all_tensor_files = list(filter(lambda file: file.endswith(TENSOR_EXT), os.listdir(tensor_folder)))
    if max_samples is not None:
        if len(all_tensor_files) < max_samples:
            logging.warning(f"{max_samples} was specified as number of samples to use but there are only "
                            f"{len(all_tensor_files)} tensor files in directory '{tensor_folder}'. Proceeding with those...")
            max_samples = len(all_tensor_files)
        tensor_files = np.random.choice(all_tensor_files, max_samples, replace=False)
    else:
        tensor_files = all_tensor_files

    num_tensor_files = len(tensor_files)
    logging.info(f"Collecting continuous stats from {num_tensor_files} of {len(all_tensor_files)} tensors at {tensor_folder}...")

    # Declare the container to hold {field_1: {sample_1: [values], sample_2: [values], field_2:...}}
    stats: DefaultDict[str, DefaultDict[str, List[float]]] = defaultdict(lambda: defaultdict(list))
    file_count = 0
    for tensor_file in tensor_files:
        _collect_continuous_stats_from_tensor_file(tensor_folder, tensor_file, stats, instances, max_arr_idx)
        file_count += 1
        if file_count % 1000 == 0:
            logging.debug(f"Collected continuous stats from {file_count}.")

    return stats, num_tensor_files


def _collect_continuous_stats_from_tensor_file(tensor_folder: str,
                                               tensor_file: str,
                                               stats: DefaultDict[str, DefaultDict[str, List[float]]],
                                               instances: List[str],
                                               max_arr_idx) -> None:
    # Inlining the method below to be able to reference more from the scope than the arguments of the function
    # 'h5py.visititems()' expects. It expects a func(<name>, <object>) => <None or return value>).
    def _field_meaning_to_values_dict(_, obj):
        if _is_continuous_valid_scalar_hd5_dataset(obj):
            value_in_tensor_file = obj[0]
            if value_in_tensor_file in CODING_VALUES_LESS_THAN_ONE:
                field_value = 0.5
            else:
                field_value = value_in_tensor_file
            dataset_name_parts = os.path.basename(obj.name).split(JOIN_CHAR)
            if len(dataset_name_parts) == 4:  # e.g. /continuous/1488_Tea-intake_0_0
                field_id = dataset_name_parts[0]
                field_meaning = dataset_name_parts[1]
                instance = dataset_name_parts[2]
                array_idx = dataset_name_parts[3]
                if instance in instances:
                    if max_arr_idx is None or (max_arr_idx is not None and int(array_idx) <= max_arr_idx):
                        stats[f"{field_meaning}{JOIN_CHAR}{field_id}{JOIN_CHAR}{instance}"][sample_id].append(field_value)
            else:  # e.g. /continuous/VentricularRate
                field_meaning = dataset_name_parts[0]
                stats[field_meaning][sample_id].append(field_value)
    tensor_file_path = os.path.join(tensor_folder, tensor_file)
    sample_id = os.path.splitext(tensor_file)[0]
    with h5py.File(tensor_file_path, 'r') as hd5_handle:
        hd5_handle.visititems(_field_meaning_to_values_dict)


def _is_continuous_valid_scalar_hd5_dataset(obj) -> bool:
    return obj.name.startswith('/continuous') and \
           isinstance(obj, h5py.Dataset) and \
           obj[0] not in CODING_VALUES_MISSING and \
           len(obj.shape) == 1
