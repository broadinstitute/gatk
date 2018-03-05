# plots.py
#
# Plotting code for Variant Filtration with Neural Nets
# This includes evaluation plots like Precision and Recall curves,
# various flavors of Reciever Operating Characteristic (ROC curves),
# As well as graphs of the metrics that are watched during neural net training.
#
# December 2016
# Sam Friedman
# sam@broadinstitute.org

# Python 2/3 friendly
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

# Imports
import os
import sys
import math
import pickle
import matplotlib
matplotlib.use('Agg') # Need this to write images from the GSA servers.  Order matters:
import matplotlib.pyplot as plt # First import matplotlib, then use Agg, then import plt

import numpy as np
from sklearn.metrics import roc_curve, auc, roc_auc_score, precision_recall_curve, average_precision_score

image_ext = '.png'

color_array = ['red', 'indigo', 'cyan', 'pink', 'purple']
key_colors = {
    'Neural Net':'green', 'CNN_SCORE':'green', 'CNN_2D':'green',
    'Heng Li Hard Filters':'lightblue',
    'GATK Hard Filters':'orange','GATK Signed Distance':'darksalmon',
    'VQSR gnomAD':'cornflowerblue', 'VQSR Single Sample':'blue', 'VQSLOD':'cornflowerblue',
    'Deep Variant':'magenta', 'QUAL':'magenta', 'DEEP_VARIANT_QUAL':'magenta',
    'Random Forest':'darkorange',
    'SNP':'cornflowerblue', 'NOT_SNP':'orange', 'INDEL':'green', 'NOT_INDEL':'red',
    'VQSLOD none':'cornflowerblue', 'VQSLOD strModel':'orange', 'VQSLOD default':'green',
    'REFERENCE':'green', 'HET_SNP':'cornflowerblue', 'HOM_SNP':'blue', 'HET_DELETION':'magenta',
    'HOM_DELETION':'violet', 'HET_INSERTION':'orange', 'HOM_INSERTION':'darkorange'
}

precision_label = 'Precision | Positive Predictive Value | TP/(TP+FP)'
recall_label = 'Recall | Sensitivity | True Positive Rate | TP/(TP+FN)'
fallout_label = 'Fallout | 1 - Specificity | False Positive Rate | FP/(FP+TN)'


def plot_precision_recall_from_scores(truth, scores, title):
    # Compute Precision-Recall and plot curve
    precision = dict()
    recall = dict()
    average_precision = dict()
    lw = 4.0
    plt.figure(figsize=(22,18))
    matplotlib.rcParams.update({'font.size': 32})
    for k in scores.keys():
        if k in key_colors:
            c = key_colors[k]
        else:
            c = np.random.choice(color_array)
        precision[k], recall[k], _ = precision_recall_curve(truth, scores[k])
        average_precision[k] = average_precision_score(truth, scores[k])
        plt.plot(recall[k], precision[k], lw=lw, color=c, label=k+' area = %0.3f' % average_precision[k])

    plt.ylim([-0.02, 1.03])
    plt.xlim([0.0, 1.0])

    plt.xlabel(recall_label)
    plt.ylabel(precision_label)
    matplotlib.rcParams.update({'font.size': 35})
    plt.title(format_title(title))

    plt.legend(loc="lower left")

    plot_name = "./figures/precision_recall_"+title+image_ext
    if not os.path.exists(os.path.dirname(plot_name)):
        os.makedirs(os.path.dirname(plot_name))
    plt.savefig(plot_name)
    print('Saved plot at:', plot_name)


def plot_simple_roc(truth, scores, title):
    fpr, tpr, _ = roc_curve(truth, scores)
    roc_score = roc_auc_score(truth, scores)
    plt.figure(figsize=(18,18))
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='Area under ROC curve = %0.2f' % roc_score)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='')
    plt.xlim([0.0, 1.0])
    plt.ylim([-0.02, 1.03])
    plt.xlabel(fallout_label)
    plt.ylabel(recall_label)
    plt.title(title.replace('_', ' ')+'\n')

    plt.legend(loc="lower right")
    plot_name = "./figures/simple_"+title+"_roc"+image_ext
    if not os.path.exists(os.path.dirname(plot_name)):
        os.makedirs(os.path.dirname(plot_name))
    plt.savefig(plot_name)
    print('Saved plot at:', plot_name)


def plot_rocs_from_scores(truth, scores, title):
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    lw = 6.0

    plt.figure(figsize=(22,18))
    matplotlib.rcParams.update({'font.size': 32})

    for k in scores.keys():
        if k in key_colors:
            c = key_colors[k]
        else:
            c = np.random.choice(color_array)

        fpr[k], tpr[k], _ = roc_curve(truth, scores[k])

        roc_auc[k] = roc_auc_score(truth, scores[k])
        plt.plot(fpr[k], tpr[k], color=c, lw=lw, label=k+' area = %0.3f' % roc_auc[k])



    plt.plot([0, 1], [0, 1], color='navy', lw=0.5, linestyle=':')
    plt.xlim([0.0, 1.0])
    plt.ylim([-0.02, 1.03])

    plt.xlabel(fallout_label)
    plt.ylabel(recall_label)

    matplotlib.rcParams.update({'font.size': 35})
    plt.title(format_title(title))

    plt.legend(loc="lower right")
    plot_name = "./figures/roc_"+title+image_ext
    if not os.path.exists(os.path.dirname(plot_name)):
        os.makedirs(os.path.dirname(plot_name))
    plt.savefig(plot_name)
    print('Saved plot at:', plot_name)


def plot_roc(model, test_data, test_truth, labels, title):
    y_pred = model.predict(test_data, batch_size=32, verbose=0)

    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i,k in enumerate(labels.keys()):
        fpr[i], tpr[i], _ = roc_curve(test_truth[:,i], y_pred[:,i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    plt.figure(figsize=(22,18))
    lw = 2
    plt.plot(fpr[0], tpr[0], color='darkorange', lw=lw, label='ROC curve (area = %0.3f)' % roc_auc[0])
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([-0.02, 1.03])
    plt.xlabel(fallout_label)
    plt.ylabel(recall_label)
    plt.title('ROC:' + str(labels) + '\n' + title)
    plt.legend(loc="lower right")
    plt.savefig("./figures/roc_"+title+image_ext)



def get_fpr_tpr_roc(model, test_data, test_truth, labels, batch_size=32, melt=False):
    y_pred = model.predict(test_data, batch_size=batch_size, verbose=0)

    if melt:
        melt_shape = (y_pred.shape[0]*y_pred.shape[1], y_pred.shape[2])
        y_pred = y_pred.reshape(melt_shape)
        test_truth = test_truth.reshape(melt_shape)

    return get_fpr_tpr_roc_pred(y_pred, test_truth, labels)


def get_fpr_tpr_roc_pred(y_pred, test_truth, labels):
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    for k in labels.keys():
        cur_idx = labels[k]
        fpr[labels[k]], tpr[labels[k]], _ = roc_curve(test_truth[:,cur_idx], y_pred[:,cur_idx])
        roc_auc[labels[k]] = auc(fpr[labels[k]], tpr[labels[k]])

    return fpr, tpr, roc_auc


def print_auc_per_class(model, test_data, test_truth, labels):
    fpr, tpr, roc_auc = get_fpr_tpr_roc(model, test_data, test_truth, labels)
    for key in labels.keys():
        print('Label %s area under ROC: %0.3f' % (key, roc_auc[labels[key]]) )


def string_auc_per_class(model, test_data, test_truth, labels):
    string = ''
    auc_sum = 0

    fpr, tpr, roc_auc = get_fpr_tpr_roc(model, test_data, test_truth, labels)
    for key in labels.keys():
        string += '\nLabel %s area under ROC: %0.3f, AUC Area Per parameter %0.8e' % (key, roc_auc[labels[key]], (roc_auc[labels[key]]/model.count_params()))
        auc_sum += roc_auc[labels[key]]

    auc_average = auc_sum/len(labels)
    aucpp_average = auc_average/model.count_params()
    string += '\n Class Averaged AUC: %0.3f, Overall AUCPP %0.8e' % (auc_average, aucpp_average)
    return string


def get_auc(model, test_data, test_truth, labels):
    mean = 0
    fpr, tpr, roc_auc = get_fpr_tpr_roc(model, test_data, test_truth, labels)
    for key in labels.keys():
        mean += roc_auc[labels[key]]

    return mean / len(labels)


def plot_roc_per_class(model, test_data, test_truth, labels, title, batch_size=32, prefix='./figures/', melt=False):
    fpr, tpr, roc_auc = get_fpr_tpr_roc(model, test_data, test_truth, labels, batch_size, melt)

    lw = 3
    plt.figure(figsize=(28,22))
    matplotlib.rcParams.update({'font.size': 34})


    for key in labels.keys():
        if key in key_colors:
            color = key_colors[key]
        else:
            color = np.random.choice(color_array)
        plt.plot( fpr[labels[key]], tpr[labels[key]], color=color, lw=lw, label=str(key)+' area under ROC: %0.3f'%roc_auc[labels[key]]  )


    plt.plot([0, 1], [0, 1], 'k:', lw=0.5)
    plt.xlim([0.0, 1.0])
    plt.ylim([-0.02, 1.03])
    plt.xlabel(fallout_label)
    plt.ylabel(recall_label)
    plt.title('ROC:'+ title + '\n')

    matplotlib.rcParams.update({'font.size': 56})
    plt.legend(loc="lower right")
    figure_path = prefix+"per_class_roc_"+title+image_ext
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    print('Saved figure at:', figure_path)


def plot_roc_per_class_predictions(predictions, test_truth, labels, title, prefix='./figures/'):
    fpr, tpr, roc_auc = get_fpr_tpr_roc_pred(predictions, test_truth, labels)
    lw = 3
    plt.figure(figsize=(28,22))
    matplotlib.rcParams.update({'font.size': 34})

    for key in labels.keys():
        if key in key_colors:
            color = key_colors[key]
        else:
            color = np.random.choice(color_array)
        plt.plot( fpr[labels[key]], tpr[labels[key]], color=color, lw=lw, label=str(key)+' area under ROC: %0.3f'%roc_auc[labels[key]]  )

    plt.plot([0, 1], [0, 1], 'k:', lw=0.5)
    plt.xlim([0.0, 1.0])
    plt.ylim([-0.02, 1.03])
    plt.xlabel(fallout_label)
    plt.ylabel(recall_label)
    plt.title('ROC:'+ title + '\n')

    matplotlib.rcParams.update({'font.size': 56})
    plt.legend(loc="lower right")
    figure_path = prefix+"per_class_roc_"+title+image_ext
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    print('Saved figure at:', figure_path)


def plot_precision_recall_per_class_predictions(predictions, truth, labels, title, prefix='./figures/'):
    # Compute Precision-Recall and plot curve
    precision = dict()
    recall = dict()
    average_precision = dict()
    lw = 4.0
    plt.figure(figsize=(22,18))
    matplotlib.rcParams.update({'font.size': 34})

    for k in labels.keys():
        if k in key_colors:
            c = key_colors[k]
        else:
            c = np.random.choice(color_array)

        precision[k], recall[k], _ = precision_recall_curve(truth[:, labels[k]], predictions[:, labels[k]])
        average_precision[k] = average_precision_score(truth[:, labels[k]], predictions[:, labels[k]])
        plt.plot(recall[k], precision[k], lw=lw, color=c, label=k+' area = %0.3f' % average_precision[k])

    plt.ylim([-0.02, 1.03])
    plt.xlim([0.0, 1.00])

    plt.xlabel(recall_label)
    plt.ylabel(precision_label)
    plt.title(title)

    plt.legend(loc="lower left")

    plot_name = prefix+"precision_recall_"+title+image_ext
    if not os.path.exists(os.path.dirname(plot_name)):
        os.makedirs(os.path.dirname(plot_name))
    plt.savefig(plot_name)
    print('Saved plot at:%s' % plot_name)


def plot_metric_history(history, title, prefix='./figures/'):
    num_plots = len ([k for k in history.history.keys() if not 'val' in k])

    row = 0
    col = 0
    rows = 4
    cols = max(2, int(math.ceil(num_plots/float(rows))))

    f, axes = plt.subplots(rows, cols, sharex=True, figsize=(36, 24))
    for k in sorted(history.history.keys()):
        if 'val' not in k:
            axes[row, col].plot(history.history[k])
            axes[row, col].set_ylabel(str(k))
            axes[row, col].set_xlabel('epoch')
            if 'val_'+k in history.history:
                axes[row, col].plot(history.history['val_'+k])
                labels = ['train', 'valid']
            else:
                labels = [k]
            axes[row, col].legend(labels, loc='upper left')

            row += 1
            if row == rows:
                row = 0
                col += 1
                if row*col >= rows*cols:
                    break

    axes[0, 1].set_title(title)
    figure_path = prefix+"metric_history_"+title+image_ext
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)


def read_tensor_to_image(args, tensor):
    eps = 1e-4
    bottom = np.amin(tensor)
    top = np.amax(tensor)

    print('Tensor bottom and top:', bottom, top , tensor.shape)

    tensor -= bottom
    tensor /= (eps+(top-bottom))
    tensor *= 254.0
    tensor = tensor.astype(np.uint8)

    if not args.channels_last:
        tensor = tensor.transpose((1,2,0))
        print('Transposed Tensor shape is:', tensor.shape)

    plt.imshow(tensor[:,:,:3])
    plt.show()


def format_title(title):
    t = title.split('/')[-1].replace('_', ' ')
    tsplit = t.split('true')
    if len(tsplit) > 1:
        t = tsplit[0]
        t += '\ntrue' + tsplit[-1]
    return t


def weight_path_to_title(wp):
    return wp.split('/')[-1].replace('__', '-').split('.')[0]

