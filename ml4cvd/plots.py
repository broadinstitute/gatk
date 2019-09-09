# plots.py

# Imports
import os
import math
import logging
import hashlib
import operator
from textwrap import wrap
from functools import reduce
from itertools import islice, product
from collections import Counter, OrderedDict, defaultdict
from typing import Iterable, DefaultDict, Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib
matplotlib.use('Agg')  # Need this to write images from the GSA servers.  Order matters:
import matplotlib.pyplot as plt  # First import matplotlib, then use Agg, then import plt
from matplotlib.ticker import NullFormatter
from matplotlib.backends.backend_pdf import PdfPages

from sklearn import manifold
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

from ml4cvd.TensorMap import TensorMap
from ml4cvd.defines import IMAGE_EXT, JOIN_CHAR, PDF_EXT

RECALL_LABEL = 'Recall | Sensitivity | True Positive Rate | TP/(TP+FN)'
FALLOUT_LABEL = 'Fallout | 1 - Specificity | False Positive Rate | FP/(FP+TN)'
PRECISION_LABEL = 'Precision | Positive Predictive Value | TP/(TP+FP)'

SUBPLOT_SIZE = 7

COLOR_ARRAY = ['tan', 'indigo', 'cyan', 'pink', 'purple', 'blue', 'chartreuse', 'deepskyblue', 'green', 'salmon', 'aqua', 'magenta', 'aquamarine', 'red',
               'coral', 'tomato', 'grey', 'black', 'maroon', 'hotpink', 'steelblue', 'orange', 'papayawhip', 'wheat', 'chocolate', 'darkkhaki', 'gold',
               'orange', 'crimson', 'slategray', 'violet', 'cadetblue', 'midnightblue', 'darkorchid', 'paleturquoise', 'plum', 'lime',
               'teal', 'peru', 'silver', 'darkgreen', 'rosybrown', 'firebrick', 'saddlebrown', 'dodgerblue', 'orangered']


def evaluate_predictions(tm: TensorMap, y_predictions: np.ndarray, y_truth: np.ndarray, title: str, folder: str, test_paths: List[str] = None,
                         max_melt: int = 5000, rocs: List[Tuple[np.ndarray, np.ndarray, Dict[str, int]]] = [],
                         scatters: List[Tuple[np.ndarray, np.ndarray, str, List[str]]] = []) -> Dict[str, float]:
    """ Evaluate predictions for a given TensorMap with truth data and plot the appropriate metrics.
    Accumulates data in the rocs and scatters lists to facilitate subplotting.

    :param tm: The TensorMap predictions to evaluate
    :param y_predictions: The predictions
    :param y_truth: The truth
    :param title: A title for the plots
    :param folder: The folder to save the plots at
    :param test_paths: The tensor paths that were predicted
    :param max_melt: For multi-dimensional prediction the maximum number of prediction to allow in the flattened array
    :param rocs: (output) List of Tuples which are inputs for ROC curve plotting to allow subplotting downstream
    :param scatters: (output) List of Tuples which are inputs for scatter plots to allow subplotting downstream
    :return: Dictionary of performance metrics with string keys for labels and float values
    """
    performance_metrics = {}
    if tm.is_categorical_any() and len(tm.shape) == 1:
        logging.info(f"For tm:{tm.name} with channel map:{tm.channel_map} examples:{y_predictions.shape[0]}")
        logging.info(f"\nSum Truth:{np.sum(y_truth, axis=0)} \nSum pred :{np.sum(y_predictions, axis=0)}")
        plot_precision_recall_per_class(y_predictions, y_truth, tm.channel_map, title, folder)
        performance_metrics.update(plot_roc_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        rocs.append((y_predictions, y_truth, tm.channel_map))
    elif tm.is_categorical() and len(tm.shape) == 2:
        melt_shape = (y_predictions.shape[0] * y_predictions.shape[1], y_predictions.shape[2])
        y_predictions = y_predictions.reshape(melt_shape)[:max_melt]
        y_truth = y_truth.reshape(melt_shape)[:max_melt]
        performance_metrics.update(plot_roc_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        performance_metrics.update(plot_precision_recall_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
    elif tm.is_categorical() and len(tm.shape) == 3:
        melt_shape = (y_predictions.shape[0] * y_predictions.shape[1] * y_predictions.shape[2], y_predictions.shape[3])
        y_predictions = y_predictions.reshape(melt_shape)[:max_melt]
        y_truth = y_truth.reshape(melt_shape)[:max_melt]
        performance_metrics.update(plot_roc_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        performance_metrics.update(plot_precision_recall_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
    elif tm.is_categorical_any() and len(tm.shape) == 4:
        melt_shape = (y_predictions.shape[0] * y_predictions.shape[1] * y_predictions.shape[2] * y_predictions.shape[3], y_predictions.shape[4])
        y_predictions = y_predictions.reshape(melt_shape)[:max_melt]
        y_truth = y_truth.reshape(melt_shape)[:max_melt]
        performance_metrics.update(plot_roc_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        performance_metrics.update(plot_precision_recall_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
    elif tm.name == 'aligned_distance':
        logging.info(f"a dist has y shape:{y_predictions.shape} and test labels has shape:{y_truth.shape}")
    elif len(tm.shape) > 1:
        prediction_flat = tm.rescale(y_predictions).flatten()
        truth_flat = tm.rescale(y_truth).flatten()
        if prediction_flat.shape[0] == truth_flat.shape[0]:
            performance_metrics.update(plot_scatter(prediction_flat, truth_flat, title, prefix=folder))
    elif tm.is_continuous():
        performance_metrics.update(plot_scatter(tm.rescale(y_predictions), tm.rescale(y_truth), title, prefix=folder, paths=test_paths))
        scatters.append((tm.rescale(y_predictions), tm.rescale(y_truth), title, test_paths))
    else:
        logging.warning(f"No evaluation clause for tensor map {tm.name}")

    if tm.name == 'median':
        plot_waves(y_predictions, y_truth, 'median_waves_' + title, folder)

    return performance_metrics


def plot_metric_history(history, title, prefix='./figures/'):
    row = 0
    col = 0
    total_plots = int(len(history.history) / 2)  # divide by 2 because we plot validation and train histories together
    cols = max(2, int(math.ceil(math.sqrt(total_plots))))
    rows = max(2, int(math.ceil(total_plots / cols)))
    f, axes = plt.subplots(rows, cols, figsize=(int(cols*SUBPLOT_SIZE), int(rows*SUBPLOT_SIZE)))
    for k in sorted(history.history.keys()):
        if 'val_' not in k:
            axes[row, col].plot(history.history[k])
            k_split = str(k).replace('output_', '').split('_')
            k_title = " ".join(OrderedDict.fromkeys(k_split))
            axes[row, col].set_title(k_title)
            axes[row, col].set_xlabel('epoch')
            if 'val_' + k in history.history:
                axes[row, col].plot(history.history['val_' + k])
                labels = ['train', 'valid']
            else:
                labels = [k]
            axes[row, col].legend(labels, loc='upper left')

            row += 1
            if row == rows:
                row = 0
                col += 1
                if col >= cols:
                    break

    plt.title(title)
    plt.tight_layout()
    figure_path = os.path.join(prefix, 'metric_history_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info(f"Saved learning curves at:{figure_path}")


def plot_scatter(prediction, truth, title, prefix='./figures/', paths=None, top_k=3, alpha=0.5):
    margin = float((np.max(truth)-np.min(truth))/100)
    plt.figure(figsize=(SUBPLOT_SIZE, SUBPLOT_SIZE))
    plt.plot([np.min(truth), np.max(truth)], [np.min(truth), np.max(truth)], linewidth=2)
    plt.plot([np.min(prediction), np.max(prediction)], [np.min(prediction), np.max(prediction)], linewidth=4)
    pearson = np.corrcoef(prediction.flatten(), truth.flatten())[1, 0]  # corrcoef returns full covariance matrix
    logging.info("Pearson coefficient is: {}".format(pearson))
    plt.scatter(prediction, truth, label=f"Pearson:{pearson:0.3f} R^2:{pearson*pearson:0.3f}", marker='.', alpha=alpha)
    if paths is not None:
        diff = np.abs(prediction-truth)
        arg_sorted = diff[:, 0].argsort()
        # The path of the best prediction, ie the inlier
        plt.text(prediction[arg_sorted[0]]+margin, truth[arg_sorted[0]]+margin, os.path.basename(paths[arg_sorted[0]]))
        # Plot the paths of the worst predictions ie the outliers
        for idx in arg_sorted[-top_k:]:
            plt.text(prediction[idx]+margin, truth[idx]+margin, os.path.basename(paths[idx]))
    plt.xlabel('Predictions')
    plt.ylabel('Actual')
    plt.title(title + '\n')
    plt.legend(loc="upper left")

    figure_path = os.path.join(prefix, 'scatter_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    logging.info("Try to save scatter plot at: {}".format(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    return {title + '_pearson': pearson}


def plot_scatters(predictions, truth, title, prefix='./figures/', paths=None, top_k=3, alpha=0.5):
    margin = float((np.max(truth) - np.min(truth)) / 100)
    plt.figure(figsize=(SUBPLOT_SIZE, SUBPLOT_SIZE))
    plt.plot([np.min(truth), np.max(truth)], [np.min(truth), np.max(truth)])
    for k in predictions:
        color = _hash_string_to_color(k)
        pearson = np.corrcoef(predictions[k].flatten(), truth.flatten())[1, 0]  # corrcoef returns full covariance matrix
        plt.plot([np.min(predictions[k]), np.max(predictions[k])], [np.min(predictions[k]), np.max(predictions[k])], color=color)
        plt.scatter(predictions[k], truth, color=color, label=str(k) + f" Pearson:{pearson:0.3f} R^2:{pearson*pearson:0.3f}", marker='.', alpha=alpha)
        if paths is not None:
            diff = np.abs(predictions[k] - truth)
            arg_sorted = diff[:, 0].argsort()
            plt.text(predictions[k][arg_sorted[0]] + margin, truth[arg_sorted[0]] + margin, os.path.basename(paths[arg_sorted[0]]))
            for idx in arg_sorted[-top_k:]:
                plt.text(predictions[k][idx] + margin, truth[idx] + margin, os.path.basename(paths[idx]))
    plt.xlabel('Predictions')
    plt.ylabel('Actual')
    plt.title(title + '\n')
    plt.legend(loc="upper left")

    figure_path = os.path.join(prefix, 'scatters_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    logging.info("Saved scatter plot at: {}".format(figure_path))


def subplot_scatters(scatters: List[Tuple[np.ndarray, np.ndarray, str, Optional[List[str]]]], prefix: str='./figures/', top_k: int=3, alpha: float=0.5):
    row = 0
    col = 0
    total_plots = len(scatters)
    cols = max(2, int(math.ceil(math.sqrt(total_plots))))
    rows = max(2, int(math.ceil(total_plots / cols)))
    fig, axes = plt.subplots(rows, cols, figsize=(cols*SUBPLOT_SIZE, rows*SUBPLOT_SIZE))
    for prediction, truth, title, paths in scatters:
        axes[row, col].plot([np.min(truth), np.max(truth)], [np.min(truth), np.max(truth)])
        axes[row, col].plot([np.min(prediction), np.max(prediction)], [np.min(prediction), np.max(prediction)])
        axes[row, col].scatter(prediction, truth, marker='.', alpha=alpha)
        margin = float((np.max(truth) - np.min(truth)) / 100)
        if paths is not None:  # If tensor paths are provided we plot the file names of top_k outliers and the #1 inlier
            diff = np.abs(prediction - truth)
            arg_sorted = diff[:, 0].argsort()
            # The path of the best prediction, ie the inlier
            axes[row, col].text(prediction[arg_sorted[0]] + margin, truth[arg_sorted[0]] + margin, os.path.basename(paths[arg_sorted[0]]))
            # Plot the paths of the worst predictions ie the outliers
            for idx in arg_sorted[-top_k:]:
                axes[row, col].text(prediction[idx] + margin, truth[idx] + margin, os.path.basename(paths[idx]))
        axes[row, col].set_xlabel('Predictions')
        axes[row, col].set_ylabel('Actual')
        axes[row, col].set_title(title + '\n')
        pearson = np.corrcoef(prediction.flatten(), truth.flatten())[1, 0]  # corrcoef returns full covariance matrix
        axes[row, col].text(0, 1, f"Pearson:{pearson:0.3f} R^2:{pearson*pearson:0.3f}", verticalalignment='bottom', transform=axes[row, col].transAxes)

        row += 1
        if row == rows:
            row = 0
            col += 1
            if col >= cols:
                break

    figure_path = prefix + 'scatters_together' + IMAGE_EXT
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    logging.info(f"Saved scatters together at: {figure_path}")


def subplot_comparison_scatters(scatters: List[Tuple[Dict[str, np.ndarray], np.ndarray, str, Optional[List[str]]]], prefix: str = './figures/', top_k: int = 3,
                                alpha: float = 0.5):
    row = 0
    col = 0
    total_plots = len(scatters)
    cols = max(2, int(math.ceil(math.sqrt(total_plots))))
    rows = max(2, int(math.ceil(total_plots / cols)))
    fig, axes = plt.subplots(rows, cols, figsize=(cols*SUBPLOT_SIZE, rows*SUBPLOT_SIZE))
    for predictions, truth, title, paths in scatters:
        for k in predictions:
            c = _hash_string_to_color(title+k)
            pearson = np.corrcoef(predictions[k].flatten(), truth.flatten())[1, 0]  # corrcoef returns full covariance matrix
            r_sqr = pearson * pearson
            axes[row, col].plot([np.min(predictions[k]), np.max(predictions[k])], [np.min(predictions[k]), np.max(predictions[k])], color=c)
            axes[row, col].scatter(predictions[k], truth, color=c, label=str(k) + ' R:%0.3f R^2:%0.3f' % (pearson, r_sqr), marker='.', alpha=alpha)
            axes[row, col].legend(loc="upper left")
            if paths is not None:  # If tensor paths are provided we plot the file names of top_k outliers and the #1 inlier
                margin = float((np.max(truth) - np.min(truth)) / 100)
                diff = np.abs(predictions[k] - truth)
                arg_sorted = diff[:, 0].argsort()
                axes[row, col].text(predictions[k][arg_sorted[0]] + margin, truth[arg_sorted[0]] + margin, os.path.basename(paths[arg_sorted[0]]))
                for idx in arg_sorted[-top_k:]:
                    axes[row, col].text(predictions[k][idx] + margin, truth[idx] + margin, os.path.basename(paths[idx]))
        axes[row, col].set_xlabel('Predictions')
        axes[row, col].set_ylabel('Actual')
        axes[row, col].set_title(title + '\n')

        row += 1
        if row == rows:
            row = 0
            col += 1
            if col >= cols:
                break

    figure_path = os.path.join(prefix, 'scatters_compared_together' + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    logging.info(f"Saved scatter comparisons together at: {figure_path}")


def plot_noise(noise):
    samples = 240
    real_weight = 2.0
    real_bias = 0.5
    x = np.linspace(10, 100, samples)
    y1_real = real_weight * x + real_bias
    y2_real = 4.0 * x + 0.8
    y1 = y1_real + (np.random.randn(*x.shape) * noise)
    y2 = y2_real + (np.random.randn(*x.shape) * noise)
    y_ratio = (y2 - y1) / y2
    y_ratio_real = (y2_real - y1_real) / y2_real
    pearson = np.corrcoef(y1.flatten(), y1_real.flatten())[1, 0]
    pearson2 = np.corrcoef(y2.flatten(), y2_real.flatten())[1, 0]
    ratio_pearson = np.corrcoef(y_ratio.flatten(), y_ratio_real.flatten())[1, 0]
    return pearson, pearson2, ratio_pearson


def plot_noisy():
    samples = 140
    p1s = []
    p2s = []
    prats = []
    noises = np.linspace(0.0, 0.01, samples)
    for n in noises:
        p1, p2, prat = plot_noise(n)
        p1s.append(1.0 - p1)
        p2s.append(1.0 - p2)
        prats.append(1.0 - prat)

    plt.figure(figsize=(28, 42))
    matplotlib.rcParams.update({'font.size': 36})
    plt.xlabel('Noise')
    plt.ylabel('Error')
    plt.scatter(noises, p1s, color='cyan', label='p1')
    plt.scatter(noises, p2s, color='green', label='p2')
    plt.scatter(noises, prats, color='red', label='p_ratio')
    plt.legend(loc="lower right")
    plt.savefig('./figures/noise_fxn.png')


def plot_value_counter(categories, counts, title, prefix='./figures/'):
    matplotlib.rcParams.update({'font.size': 14})
    counters = defaultdict(Counter)
    for k in categories:
        parts = k.split(JOIN_CHAR)
        group = parts[0]
        label = parts[1]
        counters[group][label] = counts[k]

    rows = int(math.ceil(math.sqrt(len(counters))))
    fig, axes = plt.subplots(rows, rows, figsize=(28, 24))
    for i, group in enumerate(counters):
        ax = plt.subplot(rows, rows, i + 1)
        ax.set_title(group)
        idxs = np.arange(len(counters[group]))
        ax.barh(idxs, list(counters[group].values()))
        ax.set_yticks(idxs)
        ax.set_yticklabels(list(counters[group].keys()))

    plt.tight_layout()

    figure_path = os.path.join(prefix, 'counter_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    logging.info(f"Saved counter plot at: {figure_path}")


def plot_histograms(continuous_stats, title, prefix='./figures/', num_bins=50):
    matplotlib.rcParams.update({'font.size': 14})

    rows = int(math.ceil(math.sqrt(len(continuous_stats))))
    fig, axes = plt.subplots(rows, rows, figsize=(28, 24))
    for i, group in enumerate(continuous_stats):
        a = np.array(continuous_stats[group])
        ax = plt.subplot(rows, rows, i + 1)
        ax.set_title(group + '\n Mean:%0.3f STD:%0.3f' % (np.mean(a), np.std(a)))
        ax.hist(continuous_stats[group], bins=num_bins)
    plt.tight_layout()

    figure_path = os.path.join(prefix, 'histograms_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    logging.info(f"Saved histograms plot at: {figure_path}")


def plot_histograms_in_pdf(stats: Dict[str, Dict[str, List[float]]],
                           all_samples_count: int,
                           output_file_name: str,
                           output_folder_path: str = './figures',
                           num_rows: int = 4,
                           num_cols: int = 6,
                           num_bins: int = 50,
                           title_line_width: int = 50) -> None:
    """
    Plots histograms of field values given in 'stats' in pdf
    :param stats: field names extracted from hd5 dataset names to list of values, one per sample_instance_arrayidx
    :param all_samples_count: total number of samples fields were drawn from; samples don't necessarily have values for each field
    :param output_file_name: name of output file in pdf
    :param output_folder_path: directory that output file will be written to
    :param num_rows: number of histograms that will be plotted vertically per pdf page
    :param num_cols: number of histograms that will be plotted horizontally per pdf page
    :param num_bins: number of histogram bins
    :param title_line_width: max number of characters that a plot title line will span; longer lines will be wrapped into multiple lines
    :return: None
    """
    def _sorted_chunks(d: Dict[str, Dict[str, List[float]]], size: int) -> Iterable[DefaultDict[str, List[float]]]:
        """
        :param d: dictionary to be chunked
        :param size: size of chunks
        :return: iterator of dictionary chunks with keys alphabetically sorted
        """
        it = iter(dict(sorted(d.items())))
        for i in range(0, len(d), size):
            yield {k: d[k] for k in islice(it, size)}

    subplot_width = 7.4 * num_cols
    subplot_height = 6 * num_rows
    matplotlib.rcParams.update({'font.size': 14, 'figure.figsize': (subplot_width, subplot_height)})

    figure_path = os.path.join(output_folder_path, output_file_name + PDF_EXT)
    with PdfPages(figure_path) as pdf:
        for stats_chunk in _sorted_chunks(stats, num_rows * num_cols):
            plt.subplots(num_rows, num_cols)
            for i, field in enumerate(stats_chunk):
                field_values = reduce(operator.concat, stats_chunk[field].values())
                field_sample_count = len(stats_chunk[field].keys())
                missingness = int(100 * (all_samples_count - field_sample_count) / all_samples_count)
                ax = plt.subplot(num_rows, num_cols, i + 1)
                title_text = '\n'.join(wrap(field, title_line_width))
                title_stats = '\n'.join(wrap(f"Mean:{np.mean(field_values):.2f} STD:{np.std(field_values):.2f} Missing:{missingness}% #Samples:{field_sample_count}", title_line_width))
                ax.set_title(title_text + "\n" + title_stats)
                ax.hist(field_values, bins=min(num_bins, len(set(field_values))))
            plt.tight_layout()
            pdf.savefig()

    logging.info(f"Saved histograms plot at: {figure_path}")


def plot_heatmap(stats: Dict[str, Dict[str, List[float]]],
                 output_file_name: str,
                 min_samples: int,
                 output_folder_path: str) -> None:

    """
    Plot heatmap of correlations between field pairs derived from 'stats'
    :param stats: field names extracted from hd5 dataset names to list of values, one per sample_instance_arrayidx
    :param output_file_name: name of output file in pdf
    :param output_folder_path: directory that output file will be written to
    :param min_samples: calculate correlation coefficient only if both fields have values from that many common samples
    :return: None
    """
    fields = stats.keys()
    num_fields = len(fields)
    field_pairs = product(fields, fields)
    correlations_by_field_pairs: DefaultDict[Tuple[str, str], float] = defaultdict(float)
    logging.info(f"There are {int(num_fields * (num_fields - 1) / 2)} field pairs.")
    processed_field_pair_count = 0
    nan_counter = Counter()  # keep track of if we've seen a field have NaNs
    for field1, field2 in field_pairs:
        if field1 not in nan_counter.keys() and field2 not in nan_counter.keys():
            if field1 == field2:
                correlations_by_field_pairs[(field1, field2)] = 1
            elif (field2, field1) in correlations_by_field_pairs:
                correlations_by_field_pairs[(field1, field2)] = correlations_by_field_pairs[(field2, field1)]
            else:
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
                    # Also add field2 to the counter if it has NaNs before 'continue'ing so we can skip it sooner in next iterations
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
                            correlations_by_field_pairs[(field1, field2)] = corr
                        else:
                            logging.warning(f"Pearson correlation for fields {field1} and {field2} is NaN.")
                    else:
                        logging.debug(f"Not calculating correlation for fields '{field1}' and '{field2}' "
                                      f"because they have different number of values ({len(field1_values)} vs. {len(field2_values)}).")
        else:
            continue

    logging.info(f"Total number of correlations: {len(correlations_by_field_pairs)}")

    fields_with_nans = nan_counter.keys()
    if len(fields_with_nans) != 0:
        logging.warning(f"The {len(fields_with_nans)} fields containing NaNs are: {', '.join(fields_with_nans)}.")

    # correlations_by_field_pairs = dict(random.sample(correlations_by_field_pairs.items(), 15))
    ser = pd.Series(list(correlations_by_field_pairs.values()),
                    index=pd.MultiIndex.from_tuples(correlations_by_field_pairs.keys()))
    df = ser.unstack()

    # Scale the figure size with the number of fields
    width = height = int(0.23 * num_fields)
    plt.subplots(figsize=(width, height))

    # Generate a custom diverging colourmap
    cmap = sns.diverging_palette(240, 10, n=9, as_cmap=True)

    # Plot
    ax = sns.heatmap(df, cmap=cmap, square=True, vmin=-1, vmax=1, cbar_kws={"shrink": 0.8})

    # Custom-colour missing values
    ax.set_facecolor('xkcd:light grey')

    # Scale the colourbar label font size with the number of fields
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=int(0.27 * num_fields))

    # Save figure
    fig = ax.get_figure()
    heatmap_path = os.path.join(output_folder_path, output_file_name + IMAGE_EXT)
    fig.savefig(heatmap_path)
    logging.info(f"Plotted heatmap ({df.shape[0]}x{df.shape[1]}) at: {heatmap_path}")


def plot_ecg(data, label, prefix='./figures/'):
    lw = 3
    matplotlib.rcParams.update({'font.size': 36})

    rows = int(math.ceil(math.sqrt(len(data))))
    cols = math.ceil(len(data) / rows)
    fig, axes = plt.subplots(rows, cols, figsize=(28, 24))
    for i, k in enumerate(data):
        color = _hash_string_to_color(k)
        ax = plt.subplot(rows, cols, i + 1)
        ax.set_title(k)
        ax.plot(data[k], color=color, lw=lw, label=str(k))
    plt.tight_layout()

    figure_path = os.path.join(prefix, label + '_ecg' + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    logging.info(f"Saved ECG plot at: {figure_path}")


def plot_counter(counts, title, prefix='./figures/'):
    plt.figure(figsize=(28, 32))
    matplotlib.rcParams.update({'font.size': 12})
    idxs = np.arange(len(counts))

    keyz = []
    vals = []
    for k in sorted(list(counts.keys())):
        keyz.append(k.replace('categorical/', '').replace('continuous/', ''))
        vals.append(counts[k])

    plt.barh(idxs, vals)
    plt.yticks(idxs, keyz)
    plt.tight_layout()
    plt.title(title + '\n')

    figure_path = os.path.join(prefix, 'counter_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    logging.info(f"Saved counter plot at: {figure_path}")


def plot_roc_per_class(prediction, truth, labels, title, prefix='./figures/'):
    lw = 2
    labels_to_areas = {}
    true_sums = np.sum(truth, axis=0)
    plt.figure(figsize=(SUBPLOT_SIZE*2, SUBPLOT_SIZE*2))
    fpr, tpr, roc_auc = get_fpr_tpr_roc_pred(prediction, truth, labels)

    for key in labels:
        labels_to_areas[key] = roc_auc[labels[key]]
        if 'no_' in key and len(labels) == 2:
            continue
        color = _hash_string_to_color(key)
        label_text = f'{key} area: {roc_auc[labels[key]]:.3f} n={true_sums[labels[key]]:.0f}'
        plt.plot(fpr[labels[key]], tpr[labels[key]], color=color, lw=lw, label=label_text)
        logging.info(f'ROC Label {label_text}')

    plt.xlim([0.0, 1.0])
    plt.ylim([-0.02, 1.03])
    plt.ylabel(RECALL_LABEL)
    plt.xlabel(FALLOUT_LABEL)
    plt.legend(loc="lower right")
    plt.plot([0, 1], [0, 1], 'k:', lw=0.5)
    plt.title(f'ROC {title} n={np.sum(true_sums):.0f}\n')

    figure_path = os.path.join(prefix, 'per_class_roc_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info("Saved ROC curve at: {}".format(figure_path))
    return labels_to_areas


def plot_rocs(predictions, truth, labels, title, prefix='./figures/'):
    lw = 2
    true_sums = np.sum(truth, axis=0)
    plt.figure(figsize=(SUBPLOT_SIZE*2, SUBPLOT_SIZE*2))

    for p in predictions:
        fpr, tpr, roc_auc = get_fpr_tpr_roc_pred(predictions[p], truth, labels)
        for key in labels:
            if 'no_' in key and len(labels) == 2:
                continue
            color = _hash_string_to_color(p+key)
            label_text = f'{p}_{key} area:{roc_auc[labels[key]]:.3f} n={true_sums[labels[key]]:.0f}'
            plt.plot(fpr[labels[key]], tpr[labels[key]], color=color, lw=lw, label=label_text)
            logging.info(f"ROC Label {label_text}")

    plt.xlim([0.0, 1.0])
    plt.ylim([-0.02, 1.03])
    plt.ylabel(RECALL_LABEL)
    plt.xlabel(FALLOUT_LABEL)
    plt.legend(loc='lower right')
    plt.title(f'ROC {title} n={np.sum(true_sums):.0f}\n')
    plt.plot([0, 1], [0, 1], 'k:', lw=0.5)

    figure_path = os.path.join(prefix, 'per_class_roc_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info("Saved ROC curve at: {}".format(figure_path))


def subplot_rocs(rocs: List[Tuple[np.ndarray, np.ndarray, Dict[str, int]]], prefix: str='./figures/'):
    """Log and tabulate AUCs given as nested dictionaries in the format '{model: {label: auc}}'"""
    lw = 2
    row = 0
    col = 0
    total_plots = len(rocs)
    cols = max(2, int(math.ceil(math.sqrt(total_plots))))
    rows = max(2, int(math.ceil(total_plots / cols)))
    fig, axes = plt.subplots(rows, cols, figsize=(cols*SUBPLOT_SIZE, rows*SUBPLOT_SIZE))
    for predicted, truth, labels in rocs:
        true_sums = np.sum(truth, axis=0)
        fpr, tpr, roc_auc = get_fpr_tpr_roc_pred(predicted, truth, labels)
        for key in labels:
            if 'no_' in key and len(labels) == 2:
                continue
            color = _hash_string_to_color(key)
            label_text = f'{key} area: {roc_auc[labels[key]]:.3f} n={true_sums[labels[key]]:.0f}'
            axes[row, col].plot(fpr[labels[key]], tpr[labels[key]], color=color, lw=lw, label=label_text)
            logging.info(f'ROC Label {label_text}')
        axes[row, col].set_xlim([0.0, 1.0])
        axes[row, col].set_ylim([-0.02, 1.03])
        axes[row, col].set_ylabel(RECALL_LABEL)
        axes[row, col].set_xlabel(FALLOUT_LABEL)
        axes[row, col].legend(loc='lower right')
        axes[row, col].plot([0, 1], [0, 1], 'k:', lw=0.5)
        axes[row, col].set_title(f'ROC n={np.sum(true_sums):.0f}')

        row += 1
        if row == rows:
            row = 0
            col += 1
            if col >= cols:
                break

    figure_path = prefix + 'rocs_together' + IMAGE_EXT
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)


def subplot_comparison_rocs(rocs: List[Tuple[Dict[str, np.ndarray], np.ndarray, Dict[str, int]]], prefix: str='./figures/'):
    """Log and tabulate AUCs given as nested dictionaries in the format '{model: {label: auc}}'"""
    lw = 3
    row = 0
    col = 0
    total_plots = len(rocs)
    cols = max(2, int(math.ceil(math.sqrt(total_plots))))
    rows = max(2, int(math.ceil(total_plots / cols)))
    fig, axes = plt.subplots(rows, cols, figsize=(cols*SUBPLOT_SIZE, rows*SUBPLOT_SIZE))
    for predictions, truth, labels in rocs:
        true_sums = np.sum(truth, axis=0)
        for p in predictions:
            fpr, tpr, roc_auc = get_fpr_tpr_roc_pred(predictions[p], truth, labels)
            for key in labels:
                if 'no_' in key and len(labels) == 2:
                    continue
                color = _hash_string_to_color(p + key)
                label_text = f'{p}_{key} area:{roc_auc[labels[key]]:.3f} n={true_sums[labels[key]]:.0f}'
                axes[row, col].plot(fpr[labels[key]], tpr[labels[key]], color=color, lw=lw, label=label_text)
                logging.info(f"ROC Label {label_text}")

        axes[row, col].set_xlim([0.0, 1.0])
        axes[row, col].set_ylim([-0.02, 1.03])
        axes[row, col].set_ylabel(RECALL_LABEL)
        axes[row, col].set_xlabel(FALLOUT_LABEL)
        axes[row, col].legend(loc="lower right")
        axes[row, col].plot([0, 1], [0, 1], 'k:', lw=0.5)
        axes[row, col].set_title(f'ROC n={np.sum(true_sums):.0f}\n')

        row += 1
        if row == rows:
            row = 0
            col += 1
            if col >= cols:
                break

    figure_path = os.path.join(prefix, 'rocs_compared_together' + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)


def plot_precision_recall_per_class(prediction, truth, labels, title, prefix='./figures/'):
    # Compute Precision-Recall and plot curve
    lw = 2.0
    labels_to_areas = {}
    true_sums = np.sum(truth, axis=0)
    plt.figure(figsize=(SUBPLOT_SIZE*2, SUBPLOT_SIZE*2))

    for k in labels:
        c = _hash_string_to_color(k)
        precision, recall, _ = precision_recall_curve(truth[:, labels[k]], prediction[:, labels[k]])
        average_precision = average_precision_score(truth[:, labels[k]], prediction[:, labels[k]])
        label_text = f'{k} mean precision:{average_precision:.3f} n={true_sums[labels[k]]:.0f}'
        plt.plot(recall, precision, lw=lw, color=c, label=label_text)
        logging.info(f'prAUC Label {label_text}')
        labels_to_areas[k] = average_precision

    plt.xlim([0.0, 1.00])
    plt.ylim([-0.02, 1.03])
    plt.xlabel(RECALL_LABEL)
    plt.ylabel(PRECISION_LABEL)
    plt.legend(loc="lower left")
    plt.title(f'{title} n={np.sum(true_sums):.0f}')

    figure_path = os.path.join(prefix, 'precision_recall_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info(f"Saved Precision Recall curve at: {figure_path}")
    return labels_to_areas


def plot_precision_recalls(predictions, truth, labels, title, prefix='./figures/'):
    # Compute Precision-Recall and plot curve for each model
    lw = 2.0
    true_sums = np.sum(truth, axis=0)
    plt.figure(figsize=(SUBPLOT_SIZE*2, SUBPLOT_SIZE*2))

    for p in predictions:
        for k in labels:
            c = _hash_string_to_color(p+k)
            precision, recall, _ = precision_recall_curve(truth[:, labels[k]], predictions[p][:, labels[k]])
            average_precision = average_precision_score(truth[:, labels[k]], predictions[p][:, labels[k]])
            label_text = f'{p}_{k} mean precision:{average_precision:.3f} n={true_sums[labels[k]]:.0f}'
            plt.plot(recall, precision, lw=lw, color=c, label=label_text)
            logging.info(f"prAUC Label {label_text}")

    plt.xlim([0.0, 1.00])
    plt.ylim([-0.02, 1.03])
    plt.xlabel(RECALL_LABEL)
    plt.ylabel(PRECISION_LABEL)
    plt.legend(loc="lower left")
    plt.title(f'{title} n={np.sum(true_sums):.0f}')

    figure_path = os.path.join(prefix, 'precision_recall_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info("Saved Precision Recall curve at: {}".format(figure_path))


def get_fpr_tpr_roc_pred(y_pred, test_truth, labels):
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    for k in labels.keys():
        cur_idx = labels[k]
        aser = roc_curve(test_truth[:, cur_idx], y_pred[:, cur_idx])
        fpr[labels[k]], tpr[labels[k]], _ = aser
        roc_auc[labels[k]] = auc(fpr[labels[k]], tpr[labels[k]])

    return fpr, tpr, roc_auc


def plot_waves(predicted_waves, true_waves, title, plot_path, rows=6, cols=6):
    row = 0
    col = 0
    f, axes = plt.subplots(rows, cols, sharex=True, figsize=(36, 36))
    for i in range(true_waves.shape[0]):
        axes[row, col].plot(true_waves[i, :, 0], color='blue', label='Actual Wave')
        if predicted_waves is not None:
            axes[row, col].plot(predicted_waves[i, :, 0], color='green', label='Predicted')
        axes[row, col].set_xlabel('time')
        row += 1
        if row == rows:
            row = 0
            col += 1
            if col >= cols:
                break
    plt.legend(loc="lower left")
    figure_path = os.path.join(plot_path, title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info("Saved waves at: {}".format(figure_path))


def plot_tsne(x_embed, categorical_labels, continuous_labels, gene_labels, label_dict, figure_path, alpha):
    x_embed = np.array(x_embed)
    if len(x_embed.shape) > 2:
        x_embed = np.reshape(x_embed, (x_embed.shape[0], np.prod(x_embed.shape[1:])))

    n_components = 2
    rows = max(2, len(label_dict))
    perplexities = [8, 50]
    (fig, subplots) = plt.subplots(rows, len(perplexities), figsize=(len(perplexities)*SUBPLOT_SIZE*2, rows*SUBPLOT_SIZE*2))

    p2y = {}
    for i, p in enumerate(perplexities):
        tsne = manifold.TSNE(n_components=n_components, init='pca', random_state=123, perplexity=p, learning_rate=20, n_iter_without_progress=500)
        p2y[p] = tsne.fit_transform(x_embed)

    j = -1
    for tm in label_dict:
        j += 1
        if j == rows:
            break
        categorical_subsets = {}
        if tm in categorical_labels + gene_labels:
            for c in tm.channel_map:
                categorical_subsets[tm.channel_map[c]] = label_dict[tm] == tm.channel_map[c]
        elif tm in continuous_labels:
            colors = label_dict[tm]
        for i, p in enumerate(perplexities):
            ax = subplots[j, i]
            ax.set_title(tm.name)  # + ", Perplexity=%d" % p)
            if tm in categorical_labels + gene_labels:
                color_labels = []
                for c in tm.channel_map:
                    channel_index = tm.channel_map[c]
                    color = _hash_string_to_color(c)
                    color_labels.append(c)
                    ax.scatter(p2y[p][categorical_subsets[channel_index], 0], p2y[p][categorical_subsets[channel_index], 1], c=color, alpha=alpha)
                ax.legend(color_labels, loc='lower left')
            elif tm in continuous_labels:
                points = ax.scatter(p2y[p][:, 0], p2y[p][:, 1], c=colors, alpha=alpha, cmap='jet')
                if i == len(perplexities) - 1:
                    fig.colorbar(points, ax=ax)

            ax.xaxis.set_major_formatter(NullFormatter())
            ax.yaxis.set_major_formatter(NullFormatter())
            ax.axis('tight')

    figure_path += 'tsne_plot' + IMAGE_EXT
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info(f"Saved T-SNE plot at: {figure_path}")


def _hash_string_to_color(string):
    """Hash a string to color (using hashlib and not the built-in hash for consistency between runs)"""
    return COLOR_ARRAY[int(hashlib.sha1(string.encode('utf-8')).hexdigest(), 16) % len(COLOR_ARRAY)]


if __name__ == '__main__':
    plot_noisy()
