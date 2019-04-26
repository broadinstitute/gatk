# plots.py

# Imports
import os
import math
import logging
import numpy as np
from collections import Counter, OrderedDict, defaultdict

import matplotlib
matplotlib.use('Agg')  # Need this to write images from the GSA servers.  Order matters:
import matplotlib.pyplot as plt  # First import matplotlib, then use Agg, then import plt
from sklearn.metrics import roc_curve, auc, roc_auc_score, precision_recall_curve, average_precision_score

from defines import IMAGE_EXT, JOIN_CHAR

RECALL_LABEL = 'Recall | Sensitivity | True Positive Rate | TP/(TP+FN)'
FALLOUT_LABEL = 'Fallout | 1 - Specificity | False Positive Rate | FP/(FP+TN)'
PRECISION_LABEL = 'Precision | Positive Predictive Value | TP/(TP+FP)'

COLOR_ARRAY = ['red', 'indigo', 'cyan', 'pink', 'purple', 'blue', 'chartreuse', 'darkseagreen', 'green', 'salmon',
               'magenta', 'aquamarine', 'gold', 'coral', 'tomato', 'grey', 'black', 'maroon', 'hotpink', 'steelblue',
               'orange']


def evaluate_predictions(tm, y, test_labels, test_data, plot_title, plot_folder, test_paths=None, max_melt=200000):
    performance_metrics = {}
    if tm.is_categorical_any() and len(tm.shape) == 1:
        logging.info('For tm:{} with channel map:{} examples:{}'.format(tm.name, tm.channel_map, y.shape[0]))
        logging.info('\nSum Truth:{} \nSum pred :{}'.format(np.sum(test_labels[tm.output_name()], axis=0), np.sum(y, axis=0)))
        performance_metrics.update(plot_roc_per_class(y, test_labels[tm.output_name()], tm.channel_map, plot_title, plot_folder))
    elif tm.is_categorical() and len(tm.shape) == 2:
        melt_shape = (y.shape[0]*y.shape[1], y.shape[2])
        y = y.reshape(melt_shape)[:max_melt]
        y_truth = test_labels[tm.output_name()].reshape(melt_shape)[:max_melt]
        performance_metrics.update(plot_roc_per_class(y, y_truth, tm.channel_map, plot_title, plot_folder))
        performance_metrics.update(plot_precision_recall_per_class(y, y_truth, tm.channel_map, plot_title, plot_folder))
    elif tm.is_categorical() and len(tm.shape) == 3:
        melt_shape = (y.shape[0]*y.shape[1]*y.shape[2], y.shape[3])
        y = y.reshape(melt_shape)[:max_melt]
        y_truth = test_labels[tm.output_name()].reshape(melt_shape)[:max_melt]
        performance_metrics.update(plot_roc_per_class(y, y_truth, tm.channel_map, plot_title, plot_folder))
        performance_metrics.update(plot_precision_recall_per_class(y, y_truth, tm.channel_map, plot_title, plot_folder))
    elif tm.is_categorical_any() and len(tm.shape) == 4:
        melt_shape = (y.shape[0]*y.shape[1]*y.shape[2]*y.shape[3], y.shape[4])
        y = y.reshape(melt_shape)[:max_melt]
        y_truth = test_labels[tm.output_name()].reshape(melt_shape)[:max_melt]
        performance_metrics.update(plot_roc_per_class(y, y_truth, tm.channel_map, plot_title, plot_folder))
        performance_metrics.update(plot_precision_recall_per_class(y, y_truth, tm.channel_map, plot_title, plot_folder))
    elif tm.name == 'aligned_distance':
        logging.info('a dist has y shape:{} and test labels has shape:{}'.format(y.shape, test_labels[tm.output_name()].shape))
    else:
        performance_metrics.update(plot_scatter(tm.rescale(y), tm.rescale(test_labels[tm.output_name()]), plot_title, prefix=plot_folder, paths=test_paths))

    if tm.name == 'median':
        plot_waves(y, test_labels[tm.output_name()], 'median_waves_'+plot_title, plot_folder)
        #plot_waves(None, test_data['input_strip_ecg_rest'], 'rest_waves_'+plot_title, plot_folder)

    return performance_metrics


def plot_metric_history(history, title, prefix='./figures/'):
    row = 0
    col = 0
    total_plots = int(len(history.history) / 2)  # divide by 2 because we plot validation and train histories together
    rows = max(2, int(math.ceil(math.sqrt(total_plots))))
    cols = max(2, math.ceil(total_plots / rows))
    f, axes = plt.subplots(rows, cols, sharex=True, figsize=(int(rows * 4.5), int(cols * 4.5)))

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


def plot_scatter(prediction, truth, title, prefix='./figures/', paths=None, top_k=3):
    margin = (np.max(truth)-np.min(truth))/100
    plt.figure(figsize=(16, 16))
    matplotlib.rcParams.update({'font.size': 18})
    plt.plot([np.min(truth), np.max(truth)], [np.min(truth), np.max(truth)], linewidth=2)
    plt.plot([np.min(prediction), np.max(prediction)], [np.min(prediction), np.max(prediction)], linewidth=4)
    plt.scatter(prediction, truth)
    if paths is not None:
        diff = np.abs(prediction-truth)
        argsorted = diff.argsort(axis=0)[:, 0]
        for idx in argsorted[:top_k]:
            plt.text(prediction[idx]+margin, truth[idx]+margin, os.path.basename(paths[int(idx)]))
        for idx in argsorted[-top_k:]:
            plt.text(prediction[idx]+margin, truth[idx]+margin, os.path.basename(paths[int(idx)]))
    plt.xlabel('Predictions')
    plt.ylabel('Actual')
    plt.title(title + '\n')
    pearson = np.corrcoef(prediction.flatten(), truth.flatten())[1, 0]  # corrcoef returns full covariance matrix
    logging.info("Pearson coefficient is: {}".format(pearson))
    plt.text(np.min(truth), np.max(truth), 'Pearson:%0.3f R^2:%0.3f' % (pearson, (pearson * pearson)))
    figure_path = os.path.join(prefix, 'scatter_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    logging.info("Try to save scatter plot at: {}".format(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    return {title + '_pearson': pearson}


def plot_scatters(predictions, truth, title, prefix='./figures/'):
    plt.figure(figsize=(28, 42))
    plt.rcParams.update({'font.size': 36})
    for k in predictions:
        color = COLOR_ARRAY[abs(hash(k)) % len(COLOR_ARRAY)]
        pearson = np.corrcoef(predictions[k].flatten(), truth.flatten())[1, 0]  # corrcoef returns full covariance matrix
        pearson_sqr = pearson * pearson
        plt.scatter(predictions[k], truth, color=color, label=str(k) + ' Pearson: %0.3f Pearson r^2: %0.3f' % (pearson, pearson_sqr))
    plt.xlabel('Predictions')
    plt.ylabel('Actual')
    plt.title(title + '\n')
    plt.legend(loc="lower right")

    figure_path = os.path.join(prefix, 'scatters_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    logging.info("Saved scatter plot at: {}".format(figure_path))


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


def plot_ecg(data, label, prefix='./figures/'):
    lw = 3
    matplotlib.rcParams.update({'font.size': 36})

    rows = int(math.ceil(math.sqrt(len(data))))
    cols = math.ceil(len(data) / rows)
    fig, axes = plt.subplots(rows, cols, figsize=(28, 24))
    for i, k in enumerate(data):
        color = COLOR_ARRAY[abs(hash(k)) % len(COLOR_ARRAY)]
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
    labels_to_areas = {}
    fpr, tpr, roc_auc = get_fpr_tpr_roc_pred(prediction, truth, labels)

    lw = 3
    plt.figure(figsize=(28, 22))
    matplotlib.rcParams.update({'font.size': 36})

    for key in labels:
        labels_to_areas[key] = roc_auc[labels[key]]
        if 'no_' in key and len(labels) == 2:
            continue
        color = COLOR_ARRAY[abs(hash(key)) % len(COLOR_ARRAY)]
        label_text = "{} area under ROC: {:.3f}".format(key, roc_auc[labels[key]])
        plt.plot(fpr[labels[key]], tpr[labels[key]], color=color, lw=lw, label=label_text)

    plt.plot([0, 1], [0, 1], 'k:', lw=0.5)
    plt.xlim([0.0, 1.0])
    plt.ylim([-0.02, 1.03])
    plt.xlabel(FALLOUT_LABEL)
    plt.ylabel(RECALL_LABEL)
    plt.title('ROC: ' + title + '\n')

    plt.legend(loc="lower right")
    figure_path = os.path.join(prefix, 'per_class_roc_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info("Saved ROC curve at: {}".format(figure_path))
    return labels_to_areas


def plot_rocs(predictions, truth, labels, title, prefix='./figures/'):
    lw = 3
    plt.figure(figsize=(28, 22))
    matplotlib.rcParams.update({'font.size': 36})

    for p in predictions:
        fpr, tpr, roc_auc = get_fpr_tpr_roc_pred(predictions[p], truth, labels)
        for key in labels:
            if 'no_' in key and len(labels) == 2:
                continue
            color = COLOR_ARRAY[abs(hash(p + key)) % len(COLOR_ARRAY)]
            label_text = "{}_{} area under ROC: {:.3f}".format(p, key, roc_auc[labels[key]])
            plt.plot(fpr[labels[key]], tpr[labels[key]], color=color, lw=lw, label=label_text)

    plt.plot([0, 1], [0, 1], 'k:', lw=0.5)
    plt.xlim([0.0, 1.0])
    plt.ylim([-0.02, 1.03])
    plt.xlabel(FALLOUT_LABEL)
    plt.ylabel(RECALL_LABEL)
    plt.title('ROC: ' + title + '\n')

    plt.legend(loc="lower right")
    figure_path = os.path.join(prefix, 'per_class_roc_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info("Saved ROC curve at: {}".format(figure_path))


def plot_precision_recall_per_class(prediction, truth, labels, title, prefix='./figures/'):
    # Compute Precision-Recall and plot curve
    lw = 4.0
    labels_to_areas = {}
    plt.figure(figsize=(22, 18))
    matplotlib.rcParams.update({'font.size': 34})

    for k in labels:
        c = COLOR_ARRAY[abs(hash(k)) % len(COLOR_ARRAY)]
        precision, recall, _ = precision_recall_curve(truth[:, labels[k]], prediction[:, labels[k]])
        average_precision = average_precision_score(truth[:, labels[k]], prediction[:, labels[k]])
        plt.plot(recall, precision, lw=lw, color=c, label=k + ' area = %0.3f' % average_precision)
        labels_to_areas[k] = average_precision

    plt.ylim([-0.02, 1.03])
    plt.xlim([0.0, 1.00])

    plt.xlabel(RECALL_LABEL)
    plt.ylabel(PRECISION_LABEL)
    plt.title(title)

    plt.legend(loc="lower left")
    figure_path = os.path.join(prefix, 'precision_recall_' + title + IMAGE_EXT)
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    plt.clf()
    logging.info("Saved Precision Recall curve at: {}".format(figure_path))
    return labels_to_areas


def plot_precision_recalls(predictions, truth, labels, title, prefix='./figures/'):
    # Compute Precision-Recall and plot curve for each model
    lw = 4.0
    plt.figure(figsize=(22, 18))
    matplotlib.rcParams.update({'font.size': 34})

    for p in predictions:
        for k in labels:
            c = COLOR_ARRAY[abs(hash(p + k)) % len(COLOR_ARRAY)]
            precision, recall, _ = precision_recall_curve(truth[:, labels[k]], predictions[p][:, labels[k]])
            average_precision = average_precision_score(truth[:, labels[k]], predictions[p][:, labels[k]])
            label_text = "{}_{} area under ROC: {:.3f}".format(p, k, average_precision)
            plt.plot(recall, precision, lw=lw, color=c, label=label_text)

    plt.ylim([-0.02, 1.03])
    plt.xlim([0.0, 1.00])

    plt.xlabel(RECALL_LABEL)
    plt.ylabel(PRECISION_LABEL)
    plt.title(title)

    plt.legend(loc="lower left")
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


if __name__ == '__main__':
    plot_noisy()
