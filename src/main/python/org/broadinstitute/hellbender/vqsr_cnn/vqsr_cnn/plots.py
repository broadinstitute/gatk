# plots.py
#
# Plotting code for Variant Filtration with Neural Nets
# This includes evaluation plots like Precision and Recall curves,
# various flavors of Receiver Operating Characteristic (ROC curves),
# As well as graphs of the metrics that are watched during neural net training.
#
# December 2016
# Sam Friedman
# sam@broadinstitute.org

# Imports
import os
import math
import matplotlib
import numpy as np
matplotlib.use('Agg') # Need this to write images from the GSA servers.  Order matters:
import matplotlib.pyplot as plt # First import matplotlib, then use Agg, then import plt
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


def get_fpr_tpr_roc(model, test_data, test_truth, labels, batch_size=32):
    """Get false positive and true positive rates from a classification model.

    Arguments:
        model: The model whose predictions to evaluate.
        test_data: Input testing data in the shape the model expects.
        test_truth: The true labels of the testing data
        labels: dict specifying the class labels.
        batch_size: Size of batches for prediction over the test data.

    Returns:
        dict, dict, dict: false positive rate, true positive rate, and area under ROC curve.
            The dicts all use label indices as keys. fpr and tpr dict's values are lists
            (the x and y coordinates that defines the ROC curves) and for AUC the value is a float.
    """
    y_pred = model.predict(test_data, batch_size=batch_size, verbose=0)
    return get_fpr_tpr_roc_pred(y_pred, test_truth, labels)


def get_fpr_tpr_roc_pred(y_pred, test_truth, labels):
    """Get false positive and true positive rates from predictions and true labels.

    Arguments:
        y_pred: model predictions to evaluate.
        test_truth: The true labels of the testing data
        labels: dict specifying the class labels.

    Returns:
        dict, dict, dict: false positive rate, true positive rate, and area under ROC curve.
            The dicts all use label indices as keys. fpr and tpr dict's values are lists
            (the x and y coordinates that defines the ROC curves) and for AUC the value is a float.
    """
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    for k in labels.keys():
        cur_idx = labels[k]
        fpr[labels[k]], tpr[labels[k]], _ = roc_curve(test_truth[:,cur_idx], y_pred[:,cur_idx])
        roc_auc[labels[k]] = auc(fpr[labels[k]], tpr[labels[k]])

    return fpr, tpr, roc_auc


def plot_roc_per_class(model, test_data, test_truth, labels, title, batch_size=32, prefix='./figures/'):
    """Plot a per class ROC curve.

    Arguments:
        model: The model whose predictions to evaluate.
        test_data: Input testing data in the shape the model expects.
        test_truth: The true labels of the testing data
        labels: dict specifying the class labels.
        title: the title to display on the plot.
        batch_size: Size of batches for prediction over the test data.
        prefix: path specifying where to save the plot.
    """
    fpr, tpr, roc_auc = get_fpr_tpr_roc(model, test_data, test_truth, labels, batch_size)

    lw = 3
    plt.figure(figsize=(28,22))
    matplotlib.rcParams.update({'font.size': 34})

    for key in labels.keys():
        if key in key_colors:
            color = key_colors[key]
        else:
            color = np.random.choice(color_array)
        plt.plot(fpr[labels[key]], tpr[labels[key]], color=color, lw=lw,
                 label=str(key)+' area under ROC: %0.3f'%roc_auc[labels[key]])

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


def plot_metric_history(history, title, prefix='./figures/'):
    """Plot metric history throughout training.

    Arguments:
        history: History object returned by Keras fit function.
        title: the title to display on the plot.
        prefix: path specifying where to save the plot.
    """
    num_plots = len([k for k in history.history.keys() if not 'val' in k])

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


def weight_path_to_title(wp):
    """Get a title from a model's weight path

    Arguments:
        wp: path to model's weights.

    Returns:
        str: a reformatted string
    """
    return wp.split('/')[-1].replace('__', '-').split('.')[0]

