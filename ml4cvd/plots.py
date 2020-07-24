# plots.py

# Imports
import os
import re
import math
import h5py
import glob
import logging
import hashlib
import operator
from textwrap import wrap
from functools import reduce
from datetime import datetime
from multiprocessing import Pool
from itertools import islice, product
from collections import Counter, OrderedDict, defaultdict
from typing import Iterable, DefaultDict, Dict, List, Tuple, Optional, Union, Callable

import numpy as np
import pandas as pd
from tensorflow.keras.optimizers.schedules import LearningRateSchedule

import matplotlib
matplotlib.use('Agg')  # Need this to write images from the GSA servers.  Order matters:
import matplotlib.pyplot as plt  # First import matplotlib, then use Agg, then import plt
from matplotlib.ticker import NullFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

from sklearn import manifold
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.metrics import brier_score_loss, precision_score, recall_score, f1_score
from sklearn.calibration import calibration_curve
from sksurv.metrics import concordance_index_censored
import seaborn as sns
from biosppy.signals import ecg
from scipy.ndimage.filters import gaussian_filter
from scipy import stats

from ml4cvd.TensorMap import TensorMap
from ml4cvd.metrics import concordance_index, coefficient_of_determination
from ml4cvd.defines import IMAGE_EXT, JOIN_CHAR, PDF_EXT, TENSOR_EXT, ECG_REST_LEADS, ECG_REST_MEDIAN_LEADS, PARTNERS_DATETIME_FORMAT, PARTNERS_DATE_FORMAT, HD5_GROUP_CHAR

RECALL_LABEL = 'Recall | Sensitivity | True Positive Rate | TP/(TP+FN)'
FALLOUT_LABEL = 'Fallout | 1 - Specificity | False Positive Rate | FP/(FP+TN)'
PRECISION_LABEL = 'Precision | Positive Predictive Value | TP/(TP+FP)'

SUBPLOT_SIZE = 8

COLOR_ARRAY = [
    'tan', 'indigo', 'cyan', 'pink', 'purple', 'blue', 'chartreuse', 'deepskyblue', 'green', 'salmon', 'aqua', 'magenta', 'aquamarine', 'red',
    'coral', 'tomato', 'grey', 'black', 'maroon', 'hotpink', 'steelblue', 'orange', 'papayawhip', 'wheat', 'chocolate', 'darkkhaki', 'gold',
    'orange', 'crimson', 'slategray', 'violet', 'cadetblue', 'midnightblue', 'darkorchid', 'paleturquoise', 'plum', 'lime',
    'teal', 'peru', 'silver', 'darkgreen', 'rosybrown', 'firebrick', 'saddlebrown', 'dodgerblue', 'orangered',
]

ECG_REST_PLOT_DEFAULT_YRANGE = 3.0
ECG_REST_PLOT_MAX_YRANGE = 10.0
ECG_REST_PLOT_LEADS = [
    ['strip_I','strip_aVR', 'strip_V1', 'strip_V4'],
    ['strip_II','strip_aVL', 'strip_V2', 'strip_V5'],
    ['strip_III','strip_aVF', 'strip_V3', 'strip_V6'],
]
ECG_REST_PLOT_MEDIAN_LEADS = [
    ['median_I','median_aVR', 'median_V1', 'median_V4'],
    ['median_II','median_aVL', 'median_V2', 'median_V5'],
    ['median_III','median_aVF', 'median_V3', 'median_V6'],
]
ECG_REST_PLOT_AMP_LEADS = [
    [0, 3, 6, 9],
    [1, 4, 7, 10],
    [2, 5, 8, 11],
]


def evaluate_predictions(
    tm: TensorMap, y_predictions: np.ndarray, y_truth: np.ndarray, title: str, folder: str, test_paths: List[str] = None,
    max_melt: int = 30000, rocs: List[Tuple[np.ndarray, np.ndarray, Dict[str, int]]] = [],
    scatters: List[Tuple[np.ndarray, np.ndarray, str, List[str]]] = [],
) -> Dict[str, float]:
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
    if tm.is_categorical() and tm.axes() == 1:
        logging.info(f"For tm:{tm.name} with channel map:{tm.channel_map} examples:{y_predictions.shape[0]}")
        logging.info(f"\nSum Truth:{np.sum(y_truth, axis=0)} \nSum pred :{np.sum(y_predictions, axis=0)}")
        plot_precision_recall_per_class(y_predictions, y_truth, tm.channel_map, title, folder)
        plot_prediction_calibration(y_predictions, y_truth, tm.channel_map, title, folder)
        performance_metrics.update(plot_roc_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        rocs.append((y_predictions, y_truth, tm.channel_map))
    elif tm.is_categorical() and tm.axes() == 2:
        melt_shape = (y_predictions.shape[0] * y_predictions.shape[1], y_predictions.shape[2])
        idx = np.random.choice(np.arange(melt_shape[0]), min(melt_shape[0], max_melt), replace=False)
        y_predictions = y_predictions.reshape(melt_shape)[idx]
        y_truth = y_truth.reshape(melt_shape)[idx]
        performance_metrics.update(plot_roc_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        performance_metrics.update(plot_precision_recall_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        plot_prediction_calibration(y_predictions, y_truth, tm.channel_map, title, folder)
        rocs.append((y_predictions, y_truth, tm.channel_map))
    elif tm.is_categorical() and tm.axes() == 3:
        melt_shape = (y_predictions.shape[0] * y_predictions.shape[1] * y_predictions.shape[2], y_predictions.shape[3])
        idx = np.random.choice(np.arange(melt_shape[0]), min(melt_shape[0], max_melt), replace=False)
        y_predictions = y_predictions.reshape(melt_shape)[idx]
        y_truth = y_truth.reshape(melt_shape)[idx]
        performance_metrics.update(plot_roc_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        performance_metrics.update(plot_precision_recall_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        plot_prediction_calibration(y_predictions, y_truth, tm.channel_map, title, folder)
        rocs.append((y_predictions, y_truth, tm.channel_map))
    elif tm.is_categorical() and tm.axes() == 4:
        melt_shape = (y_predictions.shape[0] * y_predictions.shape[1] * y_predictions.shape[2] * y_predictions.shape[3], y_predictions.shape[4])
        idx = np.random.choice(np.arange(melt_shape[0]), min(melt_shape[0], max_melt), replace=False)
        y_predictions = y_predictions.reshape(melt_shape)[idx]
        y_truth = y_truth.reshape(melt_shape)[idx]
        performance_metrics.update(plot_roc_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        performance_metrics.update(plot_precision_recall_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        plot_prediction_calibration(y_predictions, y_truth, tm.channel_map, title, folder)
        rocs.append((y_predictions, y_truth, tm.channel_map))
    elif tm.is_survival_curve():
        performance_metrics.update(plot_survival(y_predictions, y_truth, title, days_window=tm.days_window, prefix=folder))
        plot_survival_curves(y_predictions, y_truth, title, days_window=tm.days_window, prefix=folder, paths=test_paths)
        time_steps = tm.shape[-1]//2
        days_per_step = 1 + tm.days_window // time_steps
        predictions_at_end = 1 - np.cumprod(y_predictions[:, :time_steps], axis=-1)[:, -1]
        events_at_end = np.cumsum(y_truth[:, time_steps:], axis=-1)[:, -1]
        follow_up = np.cumsum(y_truth[:, :time_steps], axis=-1)[:, -1] * days_per_step
        logging.info(f'Shapes event {events_at_end.shape}, preds shape {predictions_at_end.shape} new ax shape {events_at_end[:, np.newaxis].shape}')
        calibration_title = f'{title}_at_{tm.days_window}_days'
        plot_prediction_calibration(predictions_at_end[:, np.newaxis], events_at_end[:, np.newaxis], {tm.name: 0}, calibration_title, folder)
        plot_survivorship(events_at_end, follow_up, predictions_at_end, tm.name, folder, tm.days_window)
    elif tm.is_time_to_event():
        c_index = concordance_index_censored(y_truth[:, 0] == 1.0, y_truth[:, 1], y_predictions[:, 0])
        concordance_return_values = ['C-Index', 'Concordant Pairs', 'Discordant Pairs', 'Tied Predicted Risk', 'Tied Event Time']
        logging.info(f"{[f'{label}: {value:.3f}' for label, value in zip(concordance_return_values, c_index)]}")
        new_title = f'{title}_C_Index_{c_index[0]:0.3f}'
        performance_metrics.update(plot_roc_per_class(y_predictions, y_truth[:, 0, np.newaxis], {f'{new_title}_vs_ROC': 0}, new_title, folder))
        calibration_title = f'{title}_at_{tm.days_window}_days'
        plot_prediction_calibration(y_predictions, y_truth[:, 0, np.newaxis], {tm.name: 0}, calibration_title, folder)
        plot_survivorship(y_truth[:, 0], y_truth[:, 1], y_predictions[:, 0], tm.name, folder, tm.days_window)
    elif tm.is_language():
        performance_metrics.update(plot_roc_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        performance_metrics.update(plot_precision_recall_per_class(y_predictions, y_truth, tm.channel_map, title, folder))
        rocs.append((y_predictions, y_truth, tm.channel_map))
    elif tm.axes() > 1 or tm.is_mesh():
        prediction_flat = tm.rescale(y_predictions).flatten()[:max_melt]
        truth_flat = tm.rescale(y_truth).flatten()[:max_melt]
        if tm.sentinel is not None:
            y_predictions = y_predictions[y_truth != tm.sentinel]
            y_truth = y_truth[y_truth != tm.sentinel]
        _plot_reconstruction(tm, y_truth, y_predictions, folder, test_paths)
        if prediction_flat.shape[0] == truth_flat.shape[0]:
            performance_metrics.update(plot_scatter(prediction_flat, truth_flat, title, prefix=folder))
    elif tm.is_continuous():
        if tm.sentinel is not None:
            y_predictions = y_predictions[y_truth != tm.sentinel, np.newaxis]
            y_truth = y_truth[y_truth != tm.sentinel, np.newaxis]
        performance_metrics.update(plot_scatter(tm.rescale(y_predictions), tm.rescale(y_truth), title, prefix=folder, paths=test_paths))
        scatters.append((tm.rescale(y_predictions), tm.rescale(y_truth), title, test_paths))
    else:
        logging.warning(f"No evaluation clause for tensor map {tm.name}")

    if tm.name == 'median':
        plot_waves(y_predictions, y_truth, 'median_waves_' + title, folder)

    return performance_metrics


def plot_metric_history(history, training_steps: int, title: str, prefix='./figures/'):
    row = 0
    col = 0
    total_plots = int(len(history.history) / 2)  # divide by 2 because we plot validation and train histories together
    cols = max(2, int(math.ceil(math.sqrt(total_plots))))
    rows = max(2, int(math.ceil(total_plots / cols)))
    f, axes = plt.subplots(rows, cols, figsize=(int(cols*SUBPLOT_SIZE), int(rows*SUBPLOT_SIZE)))
    for k in sorted(history.history.keys()):
        if not k.startswith('val_'):
            if isinstance(history.history[k][0], LearningRateSchedule):
                history.history[k] = [history.history[k][0](i * training_steps) for i in range(len(history.history[k]))]
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

    plt.tight_layout()
    figure_path = os.path.join(prefix, 'metric_history_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info(f"Saved learning curves at:{figure_path}")


def plot_rocs(predictions: Dict[str, np.ndarray], truth: np.ndarray, labels: Dict[str, int], title: str, prefix: str = './figures/'):
    """Plot Receiver Operating Characteristic (ROC) curves from a dictionary of predictions

    Typically this function is used to compare several models predictions across multiple labels.
    As a hack to avoid repetitive ROC curves for binary classification label string containing 'no_' are skipped.

    :param predictions: The keys are strings identifying the model the values are numpy arrays
                        The arrays have shape (num_samples, num_classes)
    :param truth: The true classifications of each class, one hot encoded of shape (num_samples, num_classes)
    :param labels: Dictionary mapping strings describing each class to their corresponding index in the arrays
    :param title: The name of this plot
    :param prefix: Optional path prefix where the plot will be saved
    """
    lw = 2
    true_sums = np.sum(truth, axis=0)
    plt.figure(figsize=(SUBPLOT_SIZE, SUBPLOT_SIZE))

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
    plt.plot([0, 1], [0, 1], 'k:', lw=0.5)
    plt.title(f'ROC {title} n={np.sum(true_sums):.0f}\n')

    figure_path = os.path.join(prefix, 'per_class_roc_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info(f"Saved ROC curve at: {figure_path}")


def plot_prediction_calibrations(
    predictions: Dict[str, np.ndarray], truth: np.ndarray, labels: Dict[str, int],
    title: str, prefix: str = './figures/', n_bins: int = 10,
):
    """Plot calibration performance and compute Brier Score.

    Typically this function is used to compare several models predictions across multiple labels.

    :param predictions: The keys are strings identifying the model the values are numpy arrays
                        The arrays have shape (num_samples, num_classes)
    :param truth: The true classifications of each class, one hot encoded of shape (num_samples, num_classes)
    :param labels: Dictionary mapping strings describing each class to their corresponding index in the arrays
    :param title: The name of this plot
    :param prefix: Optional path prefix where the plot will be saved
    :param n_bins: Number of bins to quantize predictions into
    """
    _ = plt.figure(figsize=(SUBPLOT_SIZE, SUBPLOT_SIZE))
    ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
    ax2 = plt.subplot2grid((3, 1), (2, 0))

    true_sums = np.sum(truth, axis=0)
    ax1.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated Brier score: 0.0")

    for p in predictions:
        for k in labels:
            color = _hash_string_to_color(p+k)
            brier_score = brier_score_loss(truth[..., labels[k]], predictions[p][..., labels[k]], pos_label=1)
            fraction_of_positives, mean_predicted_value = calibration_curve(truth[..., labels[k]], predictions[p][..., labels[k]], n_bins=n_bins)
            ax1.plot(mean_predicted_value, fraction_of_positives, "s-", label=f"{p} {k} Brier: {brier_score:0.3f}", color=color)
            ax2.hist(predictions[p][..., labels[k]], range=(0, 1), bins=10, label=f'{p} {k} n={true_sums[labels[k]]:.0f}', histtype="step", lw=2, color=color)
            logging.info(f'{p} {k} n={true_sums[labels[k]]:.0f}\nBrier score: {brier_score:0.3f}')
    ax1.set_ylabel("Fraction of positives")
    ax1.set_ylim([-0.05, 1.05])
    ax1.legend(loc="lower right")
    ax1.set_title('Calibrations plots  (reliability curve)')

    ax2.set_xlabel("Mean predicted value")
    ax2.set_ylabel("Count")
    ax2.legend(loc="upper center")
    plt.tight_layout()

    figure_path = os.path.join(prefix, 'calibrations_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    logging.info(f"Try to save calibration comparison plot at: {figure_path}")
    plt.savefig(figure_path)
    plt.clf()


def plot_prediction_calibration(
    prediction: np.ndarray, truth: np.ndarray, labels: Dict[str, int],
    title: str, prefix: str = './figures/', n_bins: int = 10,
):
    """Plot calibration performance and compute Brier Score.

    :param prediction: Array of probabilistic predictions with shape (num_samples, num_classes)
    :param truth: The true classifications of each class, one hot encoded of shape (num_samples, num_classes)
    :param labels: Dictionary mapping strings describing each class to their corresponding index in the arrays
    :param title: The name of this plot
    :param prefix: Optional path prefix where the plot will be saved
    :param n_bins: Number of bins to quantize predictions into
    """
    _, (ax1, ax3, ax2) = plt.subplots(3, figsize=(SUBPLOT_SIZE, 2*SUBPLOT_SIZE))

    true_sums = np.sum(truth, axis=0)
    ax1.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated Brier score: 0.0")
    ax3.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated Brier score: 0.0")

    for k in labels:
        y_true = truth[..., labels[k]]
        y_prob = prediction[..., labels[k]]
        color = _hash_string_to_color(k)
        brier_score = brier_score_loss(y_true, prediction[..., labels[k]], pos_label=1)
        fraction_of_positives, mean_predicted_value = calibration_curve(y_true, y_prob, n_bins=n_bins)
        ax3.plot(mean_predicted_value, fraction_of_positives, "s-", label=f"{k} Brier score: {brier_score:0.3f}", color=color)
        ax2.hist(y_prob, range=(0, 1), bins=n_bins, label=f'{k} n={true_sums[labels[k]]:.0f}', histtype="step", lw=2, color=color)

        bins = stats.mstats.mquantiles(y_prob, np.arange(0.0, 1.0, 1.0/n_bins))
        binids = np.digitize(y_prob, bins) - 1

        bin_sums = np.bincount(binids, weights=y_prob, minlength=len(bins))
        bin_true = np.bincount(binids, weights=y_true, minlength=len(bins))
        bin_total = np.bincount(binids, minlength=len(bins))

        nonzero = bin_total != 0
        prob_true = (bin_true[nonzero] / bin_total[nonzero])
        prob_pred = (bin_sums[nonzero] / bin_total[nonzero])
        ax1.plot(prob_pred, prob_true, "s-", label=f"{k} Brier score: {brier_score:0.3f}", color=color)
    ax1.set_ylabel("Fraction of positives")
    ax1.set_ylim([-0.05, 1.05])
    ax1.legend(loc="lower right")
    ax1.set_title(f'{title.replace("_", " ")}\nCalibration plot (equally sized bins)')
    ax2.set_xlabel("Mean predicted value")
    ax2.set_ylabel("Count")
    ax2.legend(loc="upper center", ncol=2)
    ax3.set_title('Calibration plot (equally spaced bins)')
    plt.tight_layout()

    figure_path = os.path.join(prefix, 'calibrations_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    logging.info(f"Try to save calibrations plot at: {figure_path}")
    plt.savefig(figure_path)
    plt.clf()


def plot_scatter(prediction, truth, title, prefix='./figures/', paths=None, top_k=3, alpha=0.5):
    margin = float((np.max(truth)-np.min(truth))/100)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(SUBPLOT_SIZE, 2 * SUBPLOT_SIZE))
    ax1.plot([np.min(truth), np.max(truth)], [np.min(truth), np.max(truth)], linewidth=2)
    ax1.plot([np.min(prediction), np.max(prediction)], [np.min(prediction), np.max(prediction)], linewidth=4)
    pearson = np.corrcoef(prediction.flatten(), truth.flatten())[1, 0]  # corrcoef returns full covariance matrix
    big_r_squared = coefficient_of_determination(truth, prediction)
    logging.info(f'Pearson:{pearson:0.3f} r^2:{pearson*pearson:0.3f} R^2:{big_r_squared:0.3f}')
    ax1.scatter(prediction, truth, label=f'Pearson:{pearson:0.3f} r^2:{pearson*pearson:0.3f} R^2:{big_r_squared:0.3f}', marker='.', alpha=alpha)
    if paths is not None:
        diff = np.abs(prediction-truth)
        arg_sorted = diff[:, 0].argsort()
        # The path of the best prediction, ie the inlier
        _text_on_plot(ax1, prediction[arg_sorted[0]]+margin, truth[arg_sorted[0]]+margin, os.path.basename(paths[arg_sorted[0]]))
        # Plot the paths of the worst predictions ie the outliers
        for idx in arg_sorted[-top_k:]:
            _text_on_plot(ax1, prediction[idx]+margin, truth[idx]+margin, os.path.basename(paths[idx]))

    ax1.set_xlabel('Predictions')
    ax1.set_ylabel('Actual')
    ax1.set_title(title + '\n')
    ax1.legend(loc="lower right")

    sns.distplot(prediction, label='Predicted', color='r', ax=ax2)
    sns.distplot(truth, label='Truth', color='b', ax=ax2)
    ax2.legend(loc="upper left")

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
        r2 = pearson*pearson
        big_r2 = coefficient_of_determination(truth.flatten(), predictions[k].flatten())
        plt.plot([np.min(predictions[k]), np.max(predictions[k])], [np.min(predictions[k]), np.max(predictions[k])], color=color)
        plt.scatter(predictions[k], truth, color=color, label=str(k) + f" Pearson:{pearson:0.3f} r^2:{r2:0.3f} R^2:{big_r2:0.3f}", marker='.', alpha=alpha)
        if paths is not None:
            diff = np.abs(predictions[k] - truth)
            arg_sorted = diff[:, 0].argsort()
            _text_on_plot(plt, predictions[k][arg_sorted[0]] + margin, truth[arg_sorted[0]] + margin, os.path.basename(paths[arg_sorted[0]]))
            for idx in arg_sorted[-top_k:]:
                _text_on_plot(plt, predictions[k][idx] + margin, truth[idx] + margin, os.path.basename(paths[idx]))
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
            _text_on_plot(axes[row, col], prediction[arg_sorted[0]] + margin, truth[arg_sorted[0]] + margin, os.path.basename(paths[arg_sorted[0]]))
            # Plot the paths of the worst predictions ie the outliers
            for idx in arg_sorted[-top_k:]:
                _text_on_plot(axes[row, col], prediction[idx] + margin, truth[idx] + margin, os.path.basename(paths[idx]))
        axes[row, col].set_xlabel('Predictions')
        axes[row, col].set_ylabel('Actual')
        axes[row, col].set_title(title + '\n')
        pearson = np.corrcoef(prediction.flatten(), truth.flatten())[1, 0]  # corrcoef returns full covariance matrix
        r2 = pearson*pearson
        big_r2 = coefficient_of_determination(truth.flatten(), prediction.flatten())
        axes[row, col].text(0, 1, f"Pearson:{pearson:0.3f} r^2:{r2:0.3f} R^2:{big_r2:0.3f}", verticalalignment='bottom', transform=axes[row, col].transAxes)

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


def subplot_comparison_scatters(
    scatters: List[Tuple[Dict[str, np.ndarray], np.ndarray, str, Optional[List[str]]]], prefix: str = './figures/', top_k: int = 3,
    alpha: float = 0.5,
):
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
            r2 = pearson * pearson
            big_r2 = coefficient_of_determination(truth.flatten(), predictions[k].flatten())
            axes[row, col].plot([np.min(predictions[k]), np.max(predictions[k])], [np.min(predictions[k]), np.max(predictions[k])], color=c)
            axes[row, col].scatter(predictions[k], truth, color=c, label=f'{k} r:{pearson:0.3f} r^2:{r2:0.3f} R^2:{big_r2:0.3f}', marker='.', alpha=alpha)
            axes[row, col].legend(loc="upper left")
            if paths is not None:  # If tensor paths are provided we plot the file names of top_k outliers and the #1 inlier
                margin = float((np.max(truth) - np.min(truth)) / 100)
                diff = np.abs(predictions[k] - truth)
                arg_sorted = diff[:, 0].argsort()
                _text_on_plot(axes[row, col], predictions[k][arg_sorted[0]] + margin, truth[arg_sorted[0]] + margin, os.path.basename(paths[arg_sorted[0]]))
                for idx in arg_sorted[-top_k:]:
                    _text_on_plot(axes[row, col], predictions[k][idx] + margin, truth[idx] + margin, os.path.basename(paths[idx]))
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


def plot_survivorship(
    events: np.ndarray, days_follow_up: np.ndarray, predictions: np.ndarray,
    title: str, prefix: str = './figures/', days_window: int = 1825,
):
    """Plot Kaplan-Meier survivorship curves and stratify by median model prediction.
    All input arrays have the same shape: (num_samples,)

    :param events: Array indicating if each sample had an event (1) or not (0) by the end of follow up
    :param days_follow_up: Array with the total days of follow up for each sample
    :param predictions: Array with model predictions of an event before the end of follow up.
    :param title: Title for the plot
    :param prefix: Path prefix where plot will be saved
    :param days_window: Maximum days of follow up
    """
    plt.figure(figsize=(SUBPLOT_SIZE, SUBPLOT_SIZE))
    days_sorted_index = np.argsort(days_follow_up)
    days_sorted = days_follow_up[days_sorted_index]
    alive_per_step = len(events)
    sick_per_step = 0
    censored = 0
    survivorship = [1.0]
    real_survivorship = [1.0]
    for cur_day, day_index in enumerate(days_sorted_index):
        if days_follow_up[day_index] > days_window:
            break
        sick_per_step += events[day_index]
        censored += 1 - events[day_index]
        alive_per_step -= events[day_index]
        survivorship.append(1 - (sick_per_step / (alive_per_step+sick_per_step)))
        real_survivorship.append(real_survivorship[cur_day] * (1 - (events[day_index] / alive_per_step)))
    logging.info(f'Cur day {cur_day} totL {len(real_survivorship)} totL {len(days_sorted)} First day {days_sorted[0]} Last day, day {days_follow_up[day_index]}, censored {censored}')
    plt.plot([0]+days_sorted[:cur_day+1], real_survivorship[:cur_day+1], marker='.', label='Survivorship')
    groups = ['High risk', 'Low risk']
    predicted_alive = {g: len(events) // 2 for g in groups}
    predicted_sick = {g: 0 for g in groups}
    predicted_days = defaultdict(list)
    predicted_survival = defaultdict(list)
    threshold = np.median(predictions)
    for cur_day, day_index in enumerate(days_sorted_index):
        if days_follow_up[day_index] > days_window:
            break
        group = 'High risk' if predictions[day_index] > threshold else 'Low risk'
        predicted_sick[group] += events[day_index]
        predicted_survival[group].append(1 - (predicted_sick[group] / (predicted_alive[group]+predicted_sick[group])))
        predicted_alive[group] -= events[day_index]
        predicted_days[group].append(days_follow_up[day_index])

    for group in groups:
        plt.plot([0]+predicted_days[group], [1]+predicted_survival[group], color='r' if 'High' in group else 'g', marker='o', label=f'{group} group had {predicted_sick[group]} events')
    plt.title(f'{title}\nEnrolled: {len(events)}, Censored: {censored:.0f}, {100 * (censored / len(events)):2.1f}%, Events: {sick_per_step:.0f}, {100 * (sick_per_step / len(events)):2.1f}%\nMax follow up: {days_window} days, {days_window // 365} years.')
    plt.xlabel('Follow up time (days)')
    plt.ylabel('Proportion Surviving')
    plt.legend(loc="lower left")

    figure_path = os.path.join(prefix, f'survivorship_fu_{days_window}_{title}{IMAGE_EXT}')
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    logging.info(f'Try to save survival plot at: {figure_path}')
    plt.savefig(figure_path)
    return {}


def plot_survival(
    prediction: np.ndarray, truth: np.ndarray, title: str, days_window: int,
    prefix: str = './figures/',
) -> Dict[str, float]:
    """Plot Kaplan-Meier survivorship and predicted proportion surviving, calculate and return C-Index

    :param prediction: Array with model predictions of an event at each time step, with shape (num_samples, intervals*2).
    :param truth: Array with survival at each time step followed by events, shape is (num_samples, intervals*2)
    :param title: Title for the plot
    :param days_window: Maximum days of follow up
    :param prefix: Path prefix where plot will be saved

    :return: Dictionary mapping metric names to their floating point values
    """
    c_index, concordant, discordant, tied_risk, tied_time = concordance_index(prediction, truth)
    logging.info(f"C-index:{c_index} concordant:{concordant} discordant:{discordant} tied_risk:{tied_risk} tied_time:{tied_time}")
    intervals = truth.shape[-1] // 2
    plt.figure(figsize=(SUBPLOT_SIZE, SUBPLOT_SIZE))

    cumulative_sick = np.cumsum(np.sum(truth[:, intervals:], axis=0))
    cumulative_censored = (truth.shape[0]-np.sum(truth[:, :intervals], axis=0))-cumulative_sick
    alive_per_step = np.sum(truth[:, :intervals], axis=0)
    sick_per_step = np.sum(truth[:, intervals:], axis=0)
    survivorship = np.cumprod(1 - (sick_per_step / alive_per_step))
    logging.info(f"Sick per step is: {sick_per_step} out of {truth.shape[0]}")
    logging.info(f"Predicted sick per step is: {list(map(int, np.sum(1-prediction[:, :intervals], axis=0)))} out of {truth.shape[0]}")
    logging.info(f"Survivors at each step is: {alive_per_step} out of {truth.shape[0]}")
    logging.info(f"Cumulative Censored: {cumulative_censored} or {np.max(truth[:, :intervals]+truth[:, intervals:])}")
    predicted_proportion = np.sum(np.cumprod(prediction[:, :intervals], axis=1), axis=0) / truth.shape[0]

    plt.plot(range(0, days_window, 1 + days_window // intervals), predicted_proportion, marker='o', label=f'Predicted Proportion C-Index:{c_index:0.3f}')
    plt.plot(range(0, days_window, 1 + days_window // intervals), survivorship, marker='o', label='Survivorship')
    plt.xlabel('Follow up time (days)')
    plt.ylabel('Proportion Surviving')
    plt.title(
        f'{title}\nEnrolled: {truth.shape[0]}, Censored: {cumulative_censored[-1]:.0f}, {100 * (cumulative_censored[-1] / truth.shape[0]):2.1f}%, '
        f'Events: {cumulative_sick[-1]:.0f}, {100 * (cumulative_sick[-1] / truth.shape[0]):2.1f}%\nMax follow up: {days_window} days, {days_window // 365} years.',
    )
    plt.legend(loc="upper right")

    figure_path = os.path.join(prefix, 'proportional_hazards_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    logging.info(f'Try to save survival plot at: {figure_path}')
    plt.savefig(figure_path)
    return {'c_index': c_index, 'concordant': concordant, 'discordant': discordant, 'tied_risk': tied_risk, 'tied_time': tied_time}


def plot_survival_curves(prediction, truth, title, days_window, prefix='./figures/', num_curves=30, paths=None):
    intervals = truth.shape[-1] // 2
    plt.figure(figsize=(SUBPLOT_SIZE*2, SUBPLOT_SIZE*2))
    predicted_survivals = np.cumprod(prediction[:, :intervals], axis=1)
    sick = np.sum(truth[:, intervals:], axis=-1)
    censor_periods = np.argmin(truth[:, :intervals], axis=-1)
    x_days = range(0, days_window, 1 + days_window // intervals)
    cur_sick = 0
    cur_healthy = 0
    min_sick = num_curves * 0.1
    for i in range(truth.shape[0]):
        p = os.path.basename(paths[i]).replace(TENSOR_EXT, "")
        if sick[i] == 1:
            sick_period = np.argmax(truth[i, intervals:])
            sick_day = sick_period*(days_window // intervals)
            plt.plot(x_days[:sick_period+2], predicted_survivals[i, :sick_period+2], label=f'Failed:{p} p:{predicted_survivals[i, sick_period]:0.2f}', color='red')
            plt.text(sick_day, predicted_survivals[i, sick_period], f'Diagnosed day:{sick_day} id:{p}')
            cur_sick += 1
            if cur_sick >= min_sick and i >= num_curves:
                break
        elif censor_periods[i] != 0:  # individual was censored before failure
            plt.plot(x_days[:censor_periods[i]], predicted_survivals[i, :censor_periods[i]], label=f'Censored:{p} p:{predicted_survivals[i, censor_periods[i]]:0.2f}', color='blue')
        elif cur_healthy < num_curves:
            plt.plot(x_days, predicted_survivals[i], label=f'Survived:{p} p:{predicted_survivals[i, -1]:0.2f}', color='green')
            cur_healthy += 1
    plt.title(title + '\n')
    plt.legend(loc="lower left")
    plt.xlabel('Follow up time (days)')
    plt.ylabel('Survival Curve Prediction')
    figure_path = os.path.join(prefix, 'survival_curves_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    logging.info("Try to save survival plot at: {}".format(figure_path))
    plt.savefig(figure_path)
    return {}


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


def plot_histograms_in_pdf(
    stats: Dict[str, Dict[str, List[float]]],
    all_samples_count: int,
    output_file_name: str,
    output_folder_path: str = './figures',
    num_rows: int = 4,
    num_cols: int = 6,
    num_bins: int = 50,
    title_line_width: int = 50,
) -> None:
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


def plot_heatmap(
    stats: Dict[str, Dict[str, List[float]]],
    output_file_name: str,
    min_samples: int,
    output_folder_path: str,
) -> None:

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
                            logging.debug(
                                f"Not calculating correlation for fields {field1} and {field2} because at least one of "
                                f"the fields has all the same values for the {num_common_samples} common samples.",
                            )
                            continue
                        corr = np.corrcoef(field1_values, field2_values)[1, 0]
                        if not math.isnan(corr):
                            correlations_by_field_pairs[(field1, field2)] = corr
                        else:
                            logging.warning(f"Pearson correlation for fields {field1} and {field2} is NaN.")
                    else:
                        logging.debug(
                            f"Not calculating correlation for fields '{field1}' and '{field2}' "
                            f"because they have different number of values ({len(field1_values)} vs. {len(field2_values)}).",
                        )
        else:
            continue

    logging.info(f"Total number of correlations: {len(correlations_by_field_pairs)}")

    fields_with_nans = nan_counter.keys()
    if len(fields_with_nans) != 0:
        logging.warning(f"The {len(fields_with_nans)} fields containing NaNs are: {', '.join(fields_with_nans)}.")

    # correlations_by_field_pairs = dict(random.sample(correlations_by_field_pairs.items(), 15))
    ser = pd.Series(
        list(correlations_by_field_pairs.values()),
        index=pd.MultiIndex.from_tuples(correlations_by_field_pairs.keys()),
    )
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


def _plot_partners_text(data: Dict[str, Union[np.ndarray, str, Dict]], fig: plt.Figure, w: float, h: float) -> None:
    # top text
    dt = datetime.strptime(data['datetime'], PARTNERS_DATETIME_FORMAT)
    dob = data['dob']
    if dob != '':
        dob = datetime.strptime(dob, PARTNERS_DATE_FORMAT)
        dob = f"{dob:%d-%b-%Y}".upper()
    age = -1
    if not np.isnan(data['age']):
        age = int(data['age'])
    sex = {value: key for key, value in data['sex'].items()}

    fig.text(0.17 / w, 8.04 / h, f"{data['lastname']}, {data['firstname']}", weight='bold')
    fig.text(3.05 / w, 8.04 / h, f"ID:{data['patientid']}", weight='bold')
    fig.text(4.56 / w, 8.04 / h, f"{dt:%d-%b-%Y %H:%M:%S}".upper(), weight='bold')
    fig.text(6.05 / w, 8.04 / h, f"{data['sitename']}", weight='bold')

    fig.text(0.17 / w, 7.77 / h, f"{dob} ({age} yr)", weight='bold')  # TODO age units
    fig.text(0.17 / w, 7.63 / h, f"{sex[1]}".title(), weight='bold')
    fig.text(0.17 / w, 7.35 / h, f"Room: ", weight='bold')  # TODO room?
    fig.text(0.17 / w, 7.21 / h, f"Loc: {data['location']}", weight='bold')

    fig.text(2.15 / w, 7.77 / h, f"Vent. rate", weight='bold')
    fig.text(2.15 / w, 7.63 / h, f"PR interval", weight='bold')
    fig.text(2.15 / w, 7.49 / h, f"QRS duration", weight='bold')
    fig.text(2.15 / w, 7.35 / h, f"QT/QTc", weight='bold')
    fig.text(2.15 / w, 7.21 / h, f"P-R-T axes", weight='bold')

    fig.text(3.91 / w, 7.77 / h, f"{int(data['rate_md'])}", weight='bold', ha='right')
    fig.text(3.91 / w, 7.63 / h, f"{int(data['pr_md'])}", weight='bold', ha='right')
    fig.text(3.91 / w, 7.49 / h, f"{int(data['qrs_md'])}", weight='bold', ha='right')
    fig.text(3.91 / w, 7.35 / h, f"{int(data['qt_md'])}/{int(data['qtc_md'])}", weight='bold', ha='right')
    fig.text(3.91 / w, 7.21 / h, f"{int(data['paxis_md'])}   {int(data['raxis_md'])}", weight='bold', ha='right')

    fig.text(4.30 / w, 7.77 / h, f"BPM", weight='bold', ha='right')
    fig.text(4.30 / w, 7.63 / h, f"ms", weight='bold', ha='right')
    fig.text(4.30 / w, 7.49 / h, f"ms", weight='bold', ha='right')
    fig.text(4.30 / w, 7.35 / h, f"ms", weight='bold', ha='right')
    fig.text(4.30 / w, 7.21 / h, f"{int(data['taxis_md'])}", weight='bold', ha='right')

    fig.text(4.75 / w, 7.21 / h, f"{data['read_md']}", wrap=True, weight='bold')

    # TODO tensorize these values from XML
    fig.text(1.28 / w, 6.65 / h, f"Technician: {''}", weight='bold')
    fig.text(1.28 / w, 6.51 / h, f"Test ind: {''}", weight='bold')
    fig.text(4.75 / w, 6.25 / h, f"Referred by: {''}", weight='bold')
    fig.text(7.63 / w, 6.25 / h, f"Electronically Signed By: {''}", weight='bold')


def _plot_partners_full(voltage: Dict[str, np.ndarray], ax: plt.Axes) -> None:
    full_voltage = np.full((12, 2500), np.nan)
    for i, lead in enumerate(voltage):
        full_voltage[i] = voltage[lead]

    # convert voltage to millivolts
    full_voltage /= 1000

    # calculate space between leads
    min_y, max_y = ax.get_ylim()
    y_offset = (max_y - min_y) / len(voltage)

    text_xoffset = 5
    text_yoffset = -0.01

    # plot signal and add labels
    for i, lead in enumerate(voltage):
        this_offset = (len(voltage) - i - 0.5) * y_offset
        ax.plot(full_voltage[i] + this_offset, color='black', linewidth=0.375)
        ax.text(
            0 + text_xoffset, this_offset + text_yoffset, lead,
            ha='left', va='top', weight='bold',
        )


def _plot_partners_clinical(voltage: Dict[str, np.ndarray], ax: plt.Axes) -> None:
    # get voltage in clinical chunks
    clinical_voltage = np.full((6, 2500), np.nan)
    halfgap = 5

    clinical_voltage[0][0:625 - halfgap] = voltage['I'][0:625 - halfgap]
    clinical_voltage[0][625 + halfgap:1250 - halfgap] = voltage['aVR'][625 + halfgap:1250 - halfgap]
    clinical_voltage[0][1250 + halfgap:1875 - halfgap] = voltage['V1'][1250 + halfgap:1875 - halfgap]
    clinical_voltage[0][1875 + halfgap:2500] = voltage['V4'][1875 + halfgap:2500]

    clinical_voltage[1][0:625 - halfgap] = voltage['II'][0:625 - halfgap]
    clinical_voltage[1][625 + halfgap:1250 - halfgap] = voltage['aVL'][625 + halfgap:1250 - halfgap]
    clinical_voltage[1][1250 + halfgap:1875 - halfgap] = voltage['V2'][1250 + halfgap:1875 - halfgap]
    clinical_voltage[1][1875 + halfgap:2500] = voltage['V5'][1875 + halfgap:2500]

    clinical_voltage[2][0:625 - halfgap] = voltage['III'][0:625 - halfgap]
    clinical_voltage[2][625 + halfgap:1250 - halfgap] = voltage['aVF'][625 + halfgap:1250 - halfgap]
    clinical_voltage[2][1250 + halfgap:1875 - halfgap] = voltage['V3'][1250 + halfgap:1875 - halfgap]
    clinical_voltage[2][1875 + halfgap:2500] = voltage['V6'][1875 + halfgap:2500]

    clinical_voltage[3] = voltage['V1']
    clinical_voltage[4] = voltage['II']
    clinical_voltage[5] = voltage['V5']

    voltage = clinical_voltage

    # convert voltage to millivolts
    voltage /= 1000

    # calculate space between leads
    min_y, max_y = ax.get_ylim()
    y_offset = (max_y - min_y) / len(voltage)

    text_xoffset = 5
    text_yoffset = -0.1

    # plot signal and add labels
    for i in range(len(voltage)):
        this_offset = (len(voltage) - i - 0.5) * y_offset
        ax.plot(voltage[i] + this_offset, color='black', linewidth=0.375)
        if i == 0:
            ax.text(
                0 + text_xoffset, this_offset + text_yoffset, 'I',
                ha='left', va='top', weight='bold',
            )
            ax.text(
                625 + text_xoffset, this_offset + text_yoffset, 'aVR',
                ha='left', va='top', weight='bold',
            )
            ax.text(
                1250 + text_xoffset, this_offset + text_yoffset, 'V1',
                ha='left', va='top', weight='bold',
            )
            ax.text(
                1875 + text_xoffset, this_offset + text_yoffset, 'V4',
                ha='left', va='top', weight='bold',
            )
        elif i == 1:
            ax.text(
                0 + text_xoffset, this_offset + text_yoffset, 'II',
                ha='left', va='top', weight='bold',
            )
            ax.text(
                625 + text_xoffset, this_offset + text_yoffset, 'aVL',
                ha='left', va='top', weight='bold',
            )
            ax.text(
                1250 + text_xoffset, this_offset + text_yoffset, 'V2',
                ha='left', va='top', weight='bold',
            )
            ax.text(
                1875 + text_xoffset, this_offset + text_yoffset, 'V5',
                ha='left', va='top', weight='bold',
            )
        elif i == 2:
            ax.text(
                0 + text_xoffset, this_offset + text_yoffset, 'III',
                ha='left', va='top', weight='bold',
            )
            ax.text(
                625 + text_xoffset, this_offset + text_yoffset, 'aVF',
                ha='left', va='top', weight='bold',
            )
            ax.text(
                1250 + text_xoffset, this_offset + text_yoffset, 'V3',
                ha='left', va='top', weight='bold',
            )
            ax.text(
                1875 + text_xoffset, this_offset + text_yoffset, 'V6',
                ha='left', va='top', weight='bold',
            )
        elif i == 3:
            ax.text(
                0 + text_xoffset, this_offset + text_yoffset, 'V1',
                ha='left', va='top', weight='bold',
            )
        elif i == 4:
            ax.text(
                0 + text_xoffset, this_offset + text_yoffset, 'II',
                ha='left', va='top', weight='bold',
            )
        elif i == 5:
            ax.text(
                0 + text_xoffset, this_offset + text_yoffset, 'V5',
                ha='left', va='top', weight='bold',
            )


def _plot_partners_figure(
        data: Dict[str, Union[np.ndarray, str, Dict]],
        plot_signal_function: Callable[[Dict[str, np.ndarray], plt.Axes], None],
        plot_mode: str,
        output_folder: str,
        run_id: str,
) -> None:
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["font.size"] = 9.5

    w, h = 11, 8.5
    fig = plt.figure(
        figsize=(w, h),
        dpi=100,
    )

    # patient info and ecg text
    _plot_partners_text(data, fig, w, h)

    # define plot area in inches
    left = 0.17
    bottom = h - 7.85
    width = w - 2 * left
    height = h - bottom - 2.3

    # ecg plot area
    ax = fig.add_axes([left/w, bottom/h, width/w, height/h])

    # voltage is in microvolts
    # the entire plot area is 5.55 inches tall, 10.66 inches wide (141 mm, 271 mm)
    # the resolution on the y-axis is 10 mm/mV
    # the resolution on the x-axis is 25 mm/s
    inch2mm = lambda inches: inches * 25.4

    # 1. set y-limit to max 14.1 mV
    y_res = 10 # mm/mV
    max_y = inch2mm(height) / y_res
    min_y = 0
    ax.set_ylim(min_y, max_y)

    # 2. set x-limit to max 10.8 s, center 10 s leads
    sampling_frequency = 250 # Hz
    x_res = 25 # mm/s
    max_x = inch2mm(width) / x_res
    x_buffer = (max_x - 10) / 2
    max_x -= x_buffer
    min_x = -x_buffer
    max_x *= sampling_frequency
    min_x *= sampling_frequency
    ax.set_xlim(min_x, max_x)

    # 3. set ticks for every 0.1 mV or every 1/25 s
    y_tick = 1 / y_res
    x_tick = 1 / x_res * sampling_frequency
    x_major_ticks = np.arange(min_x, max_x, x_tick * 5)
    x_minor_ticks = np.arange(min_x, max_x, x_tick)
    y_major_ticks = np.arange(min_y, max_y, y_tick * 5)
    y_minor_ticks = np.arange(min_y, max_y, y_tick)

    ax.set_xticks(x_major_ticks)
    ax.set_xticks(x_minor_ticks, minor=True)
    ax.set_yticks(y_major_ticks)
    ax.set_yticks(y_minor_ticks, minor=True)

    ax.tick_params(which="both", left=False, bottom=False, labelleft=False, labelbottom=False)
    ax.grid(b=True, color="r", which="major", lw=0.5)
    ax.grid(b=True, color="r", which="minor", lw=0.2)

    # signal plot
    voltage = data['2500_raw']
    plot_signal_function(voltage, ax)

    # bottom text
    fig.text(0.17 / w, 0.46 / h, f"{x_res}mm/s    {y_res}mm/mV    {sampling_frequency}Hz", ha='left', va='center', weight='bold')

    # save both pdf and png
    title = re.sub(r'[:/. ]', '', f'{plot_mode}_{data["patientid"]}_{data["datetime"]}')
    plt.savefig(os.path.join(output_folder, run_id, f'{title}{PDF_EXT}'))
    plt.savefig(os.path.join(output_folder, run_id, f'{title}{IMAGE_EXT}'))
    plt.close(fig)


def plot_partners_ecgs(args):
    plot_tensors = [
        'partners_ecg_patientid', 'partners_ecg_firstname', 'partners_ecg_lastname',
        'partners_ecg_sex',       'partners_ecg_dob',       'partners_ecg_age',
        'partners_ecg_datetime',  'partners_ecg_sitename',  'partners_ecg_location',
        'partners_ecg_read_md',   'partners_ecg_taxis_md',  'partners_ecg_rate_md',
        'partners_ecg_pr_md',     'partners_ecg_qrs_md',    'partners_ecg_qt_md',
        'partners_ecg_paxis_md',  'partners_ecg_raxis_md',  'partners_ecg_qtc_md',
    ]
    voltage_tensor = 'partners_ecg_2500_raw'
    from ml4cvd.tensor_maps_partners_ecg_labels import TMAPS
    tensor_maps_in = [TMAPS[it] for it in plot_tensors + [voltage_tensor]]
    tensor_paths = [os.path.join(args.tensors, tp) for tp in os.listdir(args.tensors) if os.path.splitext(tp)[-1].lower()==TENSOR_EXT]

    if 'clinical' == args.plot_mode:
        plot_signal_function = _plot_partners_clinical
    elif 'full' == args.plot_mode:
        plot_signal_function = _plot_partners_full
    else:
        raise ValueError(f'Unsupported plot mode: {args.plot_mode}')

    # TODO use TensorGenerator here
    # Get tensors for all hd5
    for tp in tensor_paths:
        try:
            with h5py.File(tp, 'r') as hd5:
                skip_hd5 = False
                tdict = defaultdict(dict)
                for tm in tensor_maps_in:
                    key = tm.name.split('partners_ecg_')[1]
                    try:
                        tensors = tm.tensor_from_file(tm, hd5)

                        if tm.shape[0] is not None:
                            # If not a multi-tensor tensor, wrap in array to loop through
                            tensors = np.array([tensors])
                        for i, tensor in enumerate(tensors):
                            if tm.channel_map:
                                tdict[i][key] = dict()
                                for cm in tm.channel_map:
                                    tdict[i][key][cm] = tensor[:, tm.channel_map[cm]] if tm.name == voltage_tensor else tensor[tm.channel_map[cm]]
                            else:
                                if 1 == (tm.shape[0] if tm.shape[0] is not None else tm.shape[1]):
                                    tensor = tensor.item()
                                tdict[i][key] = tensor
                    except (IndexError, KeyError, ValueError, OSError, RuntimeError) as e:
                        logging.warning(f'Could not obtain {tm.name}. Skipping plotting for all ECGs at {tp}')
                        skip_hd5 = True
                    if skip_hd5: break
                if skip_hd5: continue

                # plot each ecg
                for i in tdict:
                    _plot_partners_figure(
                        data=tdict[i],
                        plot_signal_function=plot_signal_function,
                        plot_mode=args.plot_mode,
                        output_folder=args.output_folder,
                        run_id=args.id,
                    )
        except:
            logging.exception(f"Broken tensor at: {tp}")


def plot_cross_reference(args, xref_df, title, time_description, window_start, window_end):
    # TODO make this work with multiple time windows
    if xref_df.empty:
        logging.info(f'No cross reference found for "{title}"')
        return

    title = title.replace(' ', '_')

    # compute day diffs
    day_diffs = np.array(xref_df.apply(lambda row: (row[args.time_tensor] - row[window_end]).days, axis=1))

    plt.rcParams['font.size'] = 18
    fig = plt.figure(figsize=(15,9))
    ax = fig.add_subplot(111)
    binwidth = 5
    ax.hist(day_diffs, bins=range(day_diffs.min(), day_diffs.max() + binwidth, binwidth))
    ax.set_xlabel('Days relative to event')
    ax.set_ylabel('Number of patients')
    ax.set_title(f'Distribution of {args.tensors_name} {time_description}: N={len(day_diffs)}')

    ax.text(0.05, 0.90, f'Min: {day_diffs.min()}', transform=ax.transAxes)
    ax.text(0.05, 0.85, f'Max: {day_diffs.max()}', transform=ax.transAxes)
    ax.text(0.05, 0.80, f'Median: {np.median(day_diffs):.0f}', transform=ax.transAxes)
    plt.tight_layout()

    fpath = os.path.join(args.output_folder, args.id, f'distribution_{title}{IMAGE_EXT}')
    fig.savefig(fpath)
    logging.info(f'Saved histogram of days relative to {window_end} to {fpath}')


def _ecg_rest_traces_and_text(hd5):
    """Extracts ECG resting traces from HD5 and returns a dictionary based on biosppy template"""
    path_prefix = 'ukb_ecg_rest'
    ecg_text_group = 'ecg_rest_text'
    if path_prefix not in hd5:
        raise ValueError('Tensor does not contain resting ECGs')
    twelve_leads = defaultdict(dict)
    leads = list(ECG_REST_LEADS.keys())+list(ECG_REST_MEDIAN_LEADS.keys())
    oldest_instance = sorted([instance for instance in hd5[path_prefix][leads[0]]])[0]
    for lead in leads:
        twelve_leads[lead]['raw'] = np.array(hd5[path_prefix][lead][oldest_instance], dtype=np.float32)
        if len(twelve_leads[lead]['raw']) == 5000:
            try:
                # Attempt analysis by biosppy, which may fail if not enough beats
                (
                    twelve_leads[lead]['ts_reference'], twelve_leads[lead]['filtered'], twelve_leads[lead]['rpeaks'],
                    twelve_leads[lead]['template_ts'], twelve_leads[lead]['templates'], twelve_leads[lead]['heart_rate_ts'],
                    twelve_leads[lead]['heart_rate'],
                ) = ecg.ecg(signal=twelve_leads[lead]['raw'], sampling_rate = 500., show=False)
            except:
                twelve_leads[lead]['ts_reference'] = np.linspace(0, len(twelve_leads[lead]['raw'])/500., len(twelve_leads[lead]['raw']))
    ecg_rest_text = np.array(hd5[path_prefix][ecg_text_group][oldest_instance]).tolist()
    return twelve_leads, ecg_rest_text


def _ecg_rest_ylims(yrange, yplot):
    """Returns ECG plot y-axis limits based on the range covered by ECG traces"""
    deltas   = [-1.0, 1.0]
    extremes = np.array([np.min(yplot), np.max(yplot)])
    delta_ext = extremes[1]-extremes[0]
    yrange = np.max([yrange, delta_ext*1.10])
    ylim_min = -yrange/2.0
    ylim_max = yrange/2.0
    if ((extremes[0] - ylim_min) < yrange*0.2) or \
       ((ylim_max-extremes[1]) < yrange*0.2) :
        ylim_min = extremes[0] - (yrange-delta_ext)/2.0
        ylim_max = extremes[1] + (yrange-delta_ext)/2.0
    return ylim_min, ylim_max


def _ecg_rest_yrange(twelve_leads, default_yrange, raw_scale, time_interval):
    """Returns y-range necessary not to cut any of the plotted ECG waveforms"""
    yrange = default_yrange
    for is_median, offset in zip([False, True], [3, 0]):
        for i in range(offset, offset+3):
            for j in range(0, 4):
                lead_name = ECG_REST_PLOT_LEADS[i-offset][j]
                lead = twelve_leads[lead_name]
                y_plot = np.array([elem_ * raw_scale for elem_ in lead['raw']])
                if not is_median:
                    y_plot = y_plot[
                        np.logical_and(
                            lead['ts_reference']>j*time_interval,
                            lead['ts_reference']<(j+1)*time_interval,
                        )
                    ]
                ylim_min, ylim_max = _ecg_rest_ylims(yrange, y_plot)
                yrange = ylim_max - ylim_min
    return min(yrange, ECG_REST_PLOT_MAX_YRANGE)


def _subplot_ecg_rest(twelve_leads, raw_scale, time_interval, hertz, lead_mapping, f, ax, yrange, offset, pat_df, is_median, is_blind):
    """Fills subplots with either median or raw resting ECG waveforms"""
    # plot will be in seconds vs mV, boxes are
    sec_per_box = 0.04
    mv_per_box = .1
    median_interval = 1.2  # 600 samples at 500Hz
    # if available, extract patient metadata and ECG interpretation
    if pat_df is not None:
        avl_yn = 'Y' if np.argmax(pat_df['aVL'])>0.5 else 'N'
        sl_yn  = 'Y' if np.argmax(pat_df['Sokolow_Lyon'])>0.5 else 'N'
        cor_yn = 'Y' if np.argmax(pat_df['Cornell'])>0.5 else 'N'
        sex_fm = 'F' if ((pat_df['sex'] == 'F') or (pat_df['sex'] == 'female')) else 'M'
        text   = f"ID: {pat_df['patient_id']}, sex: {sex_fm}\n"
        if not is_blind:
            text  += f"{pat_df['ecg_text']}\n"
            text  += f"LVH criteria - aVL: {avl_yn}, Sokolow-Lyon: {sl_yn}, Cornell: {cor_yn}"
        st=f.suptitle(text, x=0.0, y=1.05, ha='left', bbox=dict(facecolor='black', alpha=0.1))
    for i in range(offset, offset+3):
        for j in range(0, 4):
            lead_name = lead_mapping[i-offset][j]
            lead = twelve_leads[lead_name]
            # Convert units to mV
            if isinstance(lead, dict):
                yy = np.array([elem_ * raw_scale for elem_ in lead['raw']])
            else:
                yy = lead
            if not is_median:
                ax[i,j].set_xlim(j*time_interval,(j+1)*time_interval)
                # extract portion of waveform that is included in the actual plots
                yplot = yy[int(j*time_interval*hertz): int((j+1)*time_interval*hertz)]
            else:
                yplot = yy
            ylim_min, ylim_max = _ecg_rest_ylims(yrange, yplot)
            ax[i,j].set_ylim(ylim_min, ylim_max) # 3.0 mV range
            ax[i,j].xaxis.set_major_locator(MultipleLocator(0.2)) # major grids at every .2sec = 5 * 0.04 sec
            ax[i,j].yaxis.set_major_locator(MultipleLocator(0.5)) # major grids at every .5mV
            ax[i,j].xaxis.set_minor_locator(AutoMinorLocator(5))
            ax[i,j].yaxis.set_minor_locator(AutoMinorLocator(5))
            ax[i,j].grid(which='major', color='#CCCCCC', linestyle='--')
            ax[i,j].grid(which='minor', color='#CCCCCC', linestyle=':')
            for label in ax[i,j].xaxis.get_ticklabels()[::2]:
                label.set_visible(False)
            if len(ax[i,j].yaxis.get_ticklabels()) > 10:
                for label in ax[i,j].yaxis.get_ticklabels()[::2]:
                    label.set_visible(False)
            #normalize data in muv
            if 'ts_reference' in lead:
                ax[i,j].plot(lead['ts_reference'], yy, label='raw')
            else:
                ax[i,j].plot(np.arange(0.0, median_interval, median_interval/len(lead['raw'])), yy, label='raw')
            ax[i,j].set_title(lead_name)
            if is_median and (pat_df is not None):
                # Find where to put the R and S amp text based on ECG baseline position
                dy_ecg = (yy[-1] - ylim_min) / yrange
                if dy_ecg > 0.3: # Put in bottom right
                    dy_amp = 0.2
                else: # Put in top right
                    dy_amp = 0.85
                ax[i,j].text(0.9, dy_amp*yrange+ylim_min, f"R: {pat_df['ramp'][ECG_REST_PLOT_AMP_LEADS[i-offset][j]]:.0f}")
                ax[i,j].text(0.9, (dy_amp-0.15)*yrange+ylim_min, f"S: {pat_df['samp'][ECG_REST_PLOT_AMP_LEADS[i-offset][j]]:.0f}")


def _str_to_list_float(str_list: str) -> List[int]:
    """'[ 3. 4. nan 3 ]' --> [ 3.0, 4.0, nan, 3.0 ]"""
    tmp_str = str_list[1:-1].split()
    return list(map(float, tmp_str))


def _ecg_rest_csv_to_df(csv):
    df = pd.read_csv(csv)
    df['ramp'] = df['ramp'].apply(_str_to_list_float)
    df['samp'] = df['samp'].apply(_str_to_list_float)
    df['patient_id'] = df['patient_id'].apply(str)
    df['Sokolow_Lyon'] = df['Sokolow_Lyon'].apply(float)
    df['Cornell'] = df['Cornell'].apply(float)
    df['aVL'] = df['aVL'].apply(float)
    return df


def _remove_duplicate_rows(df, out_folder):
    arr_list = []
    pdfs = glob.glob(out_folder+'/*.pdf')
    for i, row in df.iterrows():
        if os.path.join(out_folder, row['patient_id']+'.pdf') not in pdfs:
            arr_list.append(i)
    arr = np.array(arr_list, dtype=np.int)
    return arr


def plot_ecg_rest(
    tensor_paths: List[str], rows: List[int],
    out_folder: str, is_blind: bool,
) -> None:
    """ Plots resting ECGs including annotations and LVH criteria

    :param tensor_paths: list of HDF5 file paths with ECG traces
    :param rows: indices of the subset of tensor_paths to be plotted (used by multiprocessing)
    :param out_folder: destination folder for the plots
    :param is_blind: if True, the plot gets blinded (helpful for review and annotation)
    """
    map_fields_to_tmaps = {
        'ramp': 'ecg_rest_ramplitude_raw',
        'samp': 'ecg_rest_samplitude_raw',
        'aVL': 'ecg_rest_lvh_avl',
        'Sokolow_Lyon': 'ecg_rest_lvh_sokolow_lyon',
        'Cornell': 'ecg_rest_lvh_cornell',
        }
    from ml4cvd.tensor_from_file import TMAPS
    raw_scale = 0.005 # Conversion from raw to mV
    default_yrange = ECG_REST_PLOT_DEFAULT_YRANGE # mV
    time_interval = 2.5 # time-interval per plot in seconds. ts_Reference data is in s, voltage measurement is 5 uv per lsb
    hertz = 500 # number of samples per second
    for row in rows:
        tensor_path = tensor_paths[row]
        patient_dic = {}
        patient_dic['patient_id'] = os.path.basename(tensor_path).replace(TENSOR_EXT, '')
        with h5py.File(tensor_path, 'r') as hd5:
            traces, text = _ecg_rest_traces_and_text(hd5)
            for field in map_fields_to_tmaps:
                tm = TMAPS[map_fields_to_tmaps[field]]
                patient_dic[field] = np.zeros(tm.shape)
                try:
                    patient_dic[field][:] = tm.tensor_from_file(tm, hd5)
                except ValueError as e:
                    logging.warning(e)
            is_female = 'Sex_Female_0_0' in hd5['categorical']
            patient_dic['sex'] = 'F' if is_female else 'M'
            patient_dic['ecg_text'] = text
        matplotlib.rcParams.update({'font.size': 20})
        fig, ax = plt.subplots(nrows=6, ncols=4, figsize=(24,18), tight_layout=True)
        yrange = _ecg_rest_yrange(traces, default_yrange, raw_scale, time_interval)
        _subplot_ecg_rest(
            traces, raw_scale, time_interval, hertz, ECG_REST_PLOT_LEADS, fig, ax, yrange,
            offset=3, pat_df=None, is_median=False, is_blind=is_blind,
        )
        _subplot_ecg_rest(
            traces, raw_scale, time_interval, hertz, ECG_REST_PLOT_MEDIAN_LEADS, fig, ax, yrange,
            offset=0, pat_df=patient_dic, is_median=True, is_blind=is_blind,
        )
        fig.savefig(os.path.join(out_folder, patient_dic['patient_id']+'.svg'), bbox_inches = "tight")


def plot_ecg_rest_mp(
    tensors: str,
    min_sample_id: int,
    max_sample_id: int,
    output_folder: str,
    num_workers: int = 1,
    is_blind: bool = False,
) -> None:
    """
    Generates (in parallel) plots for 12-lead resting ECGs given a CSV file pointing to the tensor HDF5s
    :param ecg_csv: name of the CSV file listing the HD5s to plot
    :param min_sample_id: smallest patient id to be plotted
    :param max_sample_id: largest patient id to be plotted
    :param output_folder: directory where output PDFs will be written to
    :param num_workers: number of parallel cores to be used
    :param is_blind: whether ECG interpretation should be included in the plot
    :return: None
    """
    tensor_paths = [os.path.join(tensors, tp) for tp in os.listdir(tensors) if os.path.splitext(tp)[-1].lower()==TENSOR_EXT]
    can_filter = True
    try:
        tensor_sample_ids = np.array([os.path.basename(tensor_path).replace(TENSOR_EXT, '') for tensor_path in tensor_paths], dtype=np.int)
        condition = (tensor_sample_ids > min_sample_id) & (tensor_sample_ids < max_sample_id)
        tensor_paths = [tensor_path for i, tensor_path in enumerate(tensor_paths) if condition[i]]
    except ValueError:
        logging.warning('Cannot select subset of tensors based on sample ids. Discarding min_ and max_sample_id')
    row_split = np.array_split(np.arange(len(tensor_paths)), num_workers)
    pool = Pool(num_workers)
    pool.starmap(plot_ecg_rest, zip([tensor_paths]*num_workers, row_split, [output_folder]*num_workers, [is_blind]*num_workers))
    pool.close()
    pool.join()


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
    plt.figure(figsize=(SUBPLOT_SIZE, SUBPLOT_SIZE))
    fpr, tpr, roc_auc = get_fpr_tpr_roc_pred(prediction, truth, labels)

    for key in labels:
        labels_to_areas[key] = roc_auc[labels[key]]
        if 'no_' in key and len(labels) == 2:
            continue
        color = _hash_string_to_color(key)
        label_text = f'{key} area: {roc_auc[labels[key]]:.3f} n={true_sums[labels[key]]:.0f}'
        plt.plot(fpr[labels[key]], tpr[labels[key]], color=color, lw=lw, label=label_text)
        logging.info(f'ROC Label {label_text} Truth shape {truth.shape}, true sums {true_sums}')

    plt.xlim([0.0, 1.0])
    plt.ylim([-0.02, 1.03])
    plt.ylabel(RECALL_LABEL)
    plt.xlabel(FALLOUT_LABEL)
    plt.legend(loc="lower right", bbox_to_anchor=(0.98, 0))
    plt.plot([0, 1], [0, 1], 'k:', lw=0.5)
    plt.title(f'ROC {title} n={truth.shape[0]:.0f}\n')

    figure_path = os.path.join(prefix, 'per_class_roc_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path, bbox_inches='tight')
    plt.clf()
    logging.info("Saved ROC curve at: {}".format(figure_path))
    return labels_to_areas


def plot_rocs(predictions, truth, labels, title, prefix='./figures/'):
    lw = 2
    true_sums = np.sum(truth, axis=0)
    plt.figure(figsize=(SUBPLOT_SIZE, SUBPLOT_SIZE))

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
    plt.plot([0, 1], [0, 1], 'k:', lw=0.5)
    plt.title(f'ROC {title} n={np.sum(true_sums):.0f}\n')

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
    plt.figure(figsize=(SUBPLOT_SIZE, SUBPLOT_SIZE))

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
    plt.figure(figsize=(SUBPLOT_SIZE, SUBPLOT_SIZE))

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

    for k in labels:
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
    perplexities = [25, 75]
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
        categorical_counts = Counter()
        if tm in categorical_labels + gene_labels:
            for c in tm.channel_map:
                categorical_subsets[tm.channel_map[c]] = label_dict[tm] == tm.channel_map[c]
                categorical_counts[tm.channel_map[c]] = np.sum(categorical_subsets[tm.channel_map[c]])
        elif tm in continuous_labels:
            colors = label_dict[tm]
        for i, p in enumerate(perplexities):
            ax = subplots[j, i]
            ax.set_title(f'{tm.name} | t-SNE perplexity:{p}')
            if tm in categorical_labels + gene_labels:
                color_labels = []
                for c in tm.channel_map:
                    channel_index = tm.channel_map[c]
                    color = _hash_string_to_color(c)
                    color_labels.append(f'{c} n={categorical_counts[tm.channel_map[c]]}')
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


def plot_find_learning_rate(
    learning_rates: List[float], losses: List[float], smoothed_losses: List[float],
    picked_learning_rate: Optional[float], figure_path: str,
):
    plt.figure(figsize=(2 * SUBPLOT_SIZE, SUBPLOT_SIZE))
    plt.title('Learning rate finder')
    cutoff = smoothed_losses[0]
    plt.ylim(min(smoothed_losses), cutoff * 1.05)
    plt.axhline(cutoff, linestyle='--', color='k', label=f'Deltas ignored above {cutoff:.2f}')
    learning_rates = np.log(learning_rates) / np.log(10)
    plt.plot(learning_rates, losses, label='Loss', c='r')
    plt.plot(learning_rates, smoothed_losses, label='Smoothed loss', c='b')
    if picked_learning_rate is not None:
        plt.axvline(np.log(picked_learning_rate) / np.log(10), label=f'Learning rate found {picked_learning_rate:.2E}', color='g', linestyle='--')
    plt.xlabel('Log_10 learning rate')
    plt.legend()
    plt.savefig(os.path.join(figure_path, f'find_learning_rate{IMAGE_EXT}'))
    plt.clf()


def plot_saliency_maps(data: np.ndarray, gradients: np.ndarray, paths: List, prefix: str):
    """Plot saliency maps of a batch of input tensors.

    Saliency maps for each input tensor in the batch will be saved at the file path indicated by prefix.
    Also creates a mean saliency map across the batch
    2D tensors are assumed to be ECGs and 3D tensors are plotted with each slice as an RGB image.
    The red channel indicates negative gradients, and the green channel positive ones.

    :param data: A batch of input tensors
    :param gradients: A corresponding batch of gradients for those inputs, must be the same shape as data
    :param paths: A List of paths corresponding to each input tensor
    :param prefix: file path prefix where saliency maps will be saved
    """
    if data.shape[-1] == 1:
        data = data[..., 0]
        gradients = gradients[..., 0]

    mean_saliency = np.zeros(data.shape[1:4] + (3,))
    for batch_i, path in enumerate(paths):
        sample_id = os.path.basename(path).replace(TENSOR_EXT, '')
        if len(data.shape) == 3:
            ecgs = {f'{sample_id}_raw': data[batch_i], 'gradients': gradients[batch_i]}
            _plot_ecgs(ecgs, f'{prefix}_{sample_id}_saliency_{batch_i}{IMAGE_EXT}')
        elif len(data.shape) == 4:
            cols = max(2, int(math.ceil(math.sqrt(data.shape[-1]))))
            rows = max(2, int(math.ceil(data.shape[-1] / cols)))
            title = f'{prefix}_{sample_id}_saliency_{batch_i}{IMAGE_EXT}'
            _plot_3d_tensor_slices_as_rgb(_saliency_map_rgb(data[batch_i], gradients[batch_i]), title, cols, rows)
            saliency = _saliency_blurred_and_scaled(gradients[batch_i], blur_radius=5.0, max_value=1.0/data.shape[0])
            mean_saliency[..., 0] -= saliency
            mean_saliency[..., 1] += saliency
        elif len(data.shape) == 5:
            for j in range(data.shape[-1]):
                cols = max(2, int(math.ceil(math.sqrt(data.shape[-2]))))
                rows = max(2, int(math.ceil(data.shape[-2] / cols)))
                name = f'{prefix}_saliency_{batch_i}_channel_{j}{IMAGE_EXT}'
                _plot_3d_tensor_slices_as_rgb(_saliency_map_rgb(data[batch_i, ..., j], gradients[batch_i, ..., j]), name, cols, rows)
                saliency = _saliency_blurred_and_scaled(gradients[batch_i, ..., j], blur_radius=5.0, max_value=1.0 / data.shape[0])
                mean_saliency[..., 0] -= saliency
                mean_saliency[..., 1] += saliency
        else:
            logging.warning(f'No method to plot saliency for data shape: {data.shape}')

    if len(data.shape) == 4:
        _plot_3d_tensor_slices_as_rgb(_scale_tensor_inplace(mean_saliency), f'{prefix}_batch_mean_saliency{IMAGE_EXT}', cols, rows)
    logging.info(f"Saved saliency maps at:{prefix}")


def _scale_tensor_inplace(tensor, min_value=0.0, max_value=1.0):
    tensor -= tensor.min()
    tensor *= (max_value - min_value) / tensor.max()
    tensor += min_value
    return tensor


def _saliency_blurred_and_scaled(gradients, blur_radius, max_value=1.0):
    blurred = gaussian_filter(gradients, sigma=blur_radius)
    _scale_tensor_inplace(blurred, max_value=max_value)
    blurred -= blurred.mean()
    return blurred


def _saliency_map_rgb(image, gradients, blur_radius=0):
    _scale_tensor_inplace(image)
    rgb_map = np.zeros(image.shape + (3,))
    blurred = _saliency_blurred_and_scaled(gradients, blur_radius)
    rgb_map[..., 0] = image - blurred
    rgb_map[..., 1] = image + blurred
    rgb_map[..., 2] = image
    rgb_map = np.clip(rgb_map, 0, 1)
    #_scale_tensor_inplace(rgb_map)
    return rgb_map


def _plot_ecgs(ecgs, figure_path, rows=3, cols=4, time_interval=2.5, raw_scale=0.005, hertz=500, lead_dictionary=ECG_REST_LEADS):
    index2leads = {v: k for k, v in lead_dictionary.items()}
    _, axes = plt.subplots(rows, cols, figsize=(18, 16), sharey=True)
    for i in range(rows):
        for j in range(cols):
            start = int(i*time_interval*hertz)
            stop = int((i+1)*time_interval*hertz)
            axes[i, j].set_xlim(start, stop)
            for label in ecgs:
                axes[i, j].plot(range(start, stop), ecgs[label][start:stop, j + i*cols] * raw_scale, label=label)
            axes[i, j].legend(loc='lower right')
            axes[i, j].set_xlabel('milliseconds')
            axes[i, j].set_ylabel('mV')
            axes[i, j].set_title(index2leads[j + i*cols])
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()


def _plot_3d_tensor_slices_as_rgb(tensor, figure_path, cols=3, rows=10):
    _, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4))
    for i in range(tensor.shape[-2]):
        axes[i // cols, i % cols].imshow(tensor[:, :, i, :])
        axes[i // cols, i % cols].set_yticklabels([])
        axes[i // cols, i % cols].set_xticklabels([])

    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()


def _hash_string_to_color(string):
    """Hash a string to color (using hashlib and not the built-in hash for consistency between runs)"""
    return COLOR_ARRAY[int(hashlib.sha1(string.encode('utf-8')).hexdigest(), 16) % len(COLOR_ARRAY)]


def _text_on_plot(axes, x, y, text, alpha=0.8, background='white'):
    t = axes.text(x, y, text)
    t.set_bbox({'facecolor': background, 'alpha': alpha, 'edgecolor': background})


def _plot_reconstruction(
        tm: TensorMap, y_true: np.ndarray, y_pred: np.ndarray,
        folder: str, paths: List[str],
):
    num_samples = 3
    logging.info(f'Plotting {num_samples} reconstructions of {tm}.')
    if None in tm.shape:  # can't handle dynamic shapes
        return
    for i in range(num_samples):
        title = f'{tm.name}_{os.path.basename(paths[i]).replace(TENSOR_EXT, "")}_reconstruction'
        y = y_true[i].reshape(tm.shape)
        yp = y_pred[i].reshape(tm.shape)
        if tm.axes() == 2:
            fig = plt.figure(figsize=(SUBPLOT_SIZE, SUBPLOT_SIZE * num_samples))
            for j in range(tm.shape[1]):
                plt.subplot(tm.shape[1], 1, j + 1)
                plt.plot(y[:, j], c='k', linestyle='--', label='original')
                plt.plot(yp[:, j], c='b', label='reconstruction')
                if j == 0:
                    plt.title(title)
                    plt.legend()
            plt.tight_layout()
        # TODO: implement 3d, 4d
        plt.savefig(os.path.join(folder, title + IMAGE_EXT))
        plt.clf()


if __name__ == '__main__':
    plot_noisy()
