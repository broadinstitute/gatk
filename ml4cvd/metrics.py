# metrics.py
import logging
import numpy as np
import tensorflow as tf
import keras.backend as K

from sklearn.metrics import roc_curve, auc, average_precision_score

from keras.losses import binary_crossentropy, categorical_crossentropy, logcosh, cosine_proximity, mean_squared_error, mean_absolute_error

STRING_METRICS = ['categorical_crossentropy','binary_crossentropy','mean_absolute_error','mae',
                  'mean_squared_error', 'mse', 'cosine_proximity', 'logcosh']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~ Metrics ~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def weighted_crossentropy(weights, name='anonymous'):
    """A weighted version of keras.objectives.categorical_crossentropy
    
    Arguments:
        weights = np.array([0.5,2,10]) # Class one at 0.5, class 2 twice the normal weights, class 3 10x.
        name: string identifying the loss to differentiate when models have multiple losses 
    
    Returns:
        keras loss function named name+'_weighted_loss'
    
    """
    string_globe = 'global ' + name + '_weights\n'
    string_globe += 'global ' + name + '_kweights\n'
    string_globe += name + '_weights = weights\n'
    string_globe += name + '_kweights = K.variable('+name+'_weights)\n'
    exec(string_globe, globals(), locals())
    fxn_postfix = '_weighted_loss'
    string_fxn = 'def '+ name + fxn_postfix + '(y_true, y_pred):\n'
    string_fxn += '\ty_pred /= K.sum(y_pred, axis=-1, keepdims=True)\n'
    string_fxn += '\ty_pred = K.clip(y_pred, K.epsilon(), 1 - K.epsilon())\n'
    string_fxn += '\tloss = y_true * K.log(y_pred) * ' + name + '_kweights\n'
    string_fxn += '\tloss = -K.sum(loss, -1)\n'
    string_fxn += '\treturn loss\n'
    exec(string_fxn, globals(), locals())
    loss_fxn = eval(name + fxn_postfix, globals(), locals())
    return loss_fxn


def angle_between_batches(tensors):
    l0 = K.sqrt(K.sum(K.square(tensors[0]), axis=-1, keepdims=True) + K.epsilon())
    l1 = K.sqrt(K.sum(K.square(tensors[1]), axis=-1, keepdims=True) + K.epsilon())
    numerator = K.sum(tensors[0]*tensors[1], axis=-1, keepdims=True)
    return tf.acos(numerator / (l0*l1))


def sum_pred_loss(y_true, y_pred):
    return K.sum(y_pred, axis=-1)


def two_batch_euclidean(tensors):
    return K.sqrt(K.sum(K.square(tensors[0] - tensors[1]), axis=-1, keepdims=True) + K.epsilon())


def custom_loss_keras(user_id, encodings):
    pairwise_diff = K.expand_dims(encodings, 0) - K.expand_dims(encodings, 1)
    pairwise_squared_distance = K.sum(K.square(pairwise_diff), axis=-1)
    pairwise_distance = K.sqrt(pairwise_squared_distance + K.epsilon())

    user_id = K.squeeze(user_id, axis=1)  # remove the axis added by Keras
    pairwise_equal = K.equal(K.expand_dims(user_id, 0), K.expand_dims(user_id, 1))

    pos_neg = K.cast(pairwise_equal, K.floatx()) * 2 - 1
    return K.sum(pairwise_distance * pos_neg, axis=-1) / 2


def euclid_dist(v):
    return (v[0] - v[1])**2

def angle_between_batches(tensors):
    l0 = K.sqrt(K.sum(K.square(tensors[0]), axis=-1, keepdims=True) + K.epsilon())
    l1 = K.sqrt(K.sum(K.square(tensors[1]), axis=-1, keepdims=True) + K.epsilon())
    numerator = K.sum(tensors[0]*tensors[1], axis=-1, keepdims=True)
    return tf.acos(numerator / (l0*l1))


def paired_angle_between_batches(tensors):
    l0 = K.sqrt(K.sum(K.square(tensors[0]), axis=-1, keepdims=True) + K.epsilon())
    l1 = K.sqrt(K.sum(K.square(tensors[1]), axis=-1, keepdims=True) + K.epsilon())
    numerator = K.sum(tensors[0]*tensors[1], axis=-1, keepdims=True)
    angle_w_self = tf.acos(numerator / (l0*l1))
    # This is very hacky! we assume batch sizes are odd and reverse the batch to compare to others.
    l1_other = K.sqrt(K.sum(K.square(K.reverse(tensors[1], 0)), axis=-1, keepdims=True) + K.epsilon())
    other_numerator = K.sum(tensors[0]*K.reverse(tensors[1], 0), axis=-1, keepdims=True)
    angle_w_other = tf.acos(other_numerator / (l0*l1_other))
    return angle_w_self - angle_w_other


def ignore_zeros_l2(y_true, y_pred):
    mask = K.cast(K.not_equal(y_true, 0), K.floatx())
    return mean_squared_error(y_true * mask, y_pred * mask)


def ignore_zeros_logcosh(y_true, y_pred):
    mask = K.cast(K.not_equal(y_true, 0), K.floatx())
    return logcosh(y_true * mask, y_pred * mask)


def sentinel_logcosh_loss(sentinel: float):
    def ignore_sentinel_logcosh(y_true, y_pred):
        mask = K.cast(K.not_equal(y_true, sentinel), K.floatx())
        return logcosh(y_true * mask, y_pred * mask)
    return ignore_sentinel_logcosh


def sum_pred_loss(y_true, y_pred):
    return K.sum(y_pred, axis=-1)


def two_batch_euclidean(tensors):
    return K.sqrt(K.sum(K.square(tensors[0] - tensors[1]), axis=-1, keepdims=True) + K.epsilon())


def custom_loss_keras(user_id, encodings):
    pairwise_diff = K.expand_dims(encodings, 0) - K.expand_dims(encodings, 1)
    pairwise_squared_distance = K.sum(K.square(pairwise_diff), axis=-1)
    pairwise_distance = K.sqrt(pairwise_squared_distance + K.epsilon())

    user_id = K.squeeze(user_id, axis=1)  # remove the axis added by Keras
    pairwise_equal = K.equal(K.expand_dims(user_id, 0), K.expand_dims(user_id, 1))

    pos_neg = K.cast(pairwise_equal, K.floatx()) * 2 - 1
    return K.sum(pairwise_distance * pos_neg, axis=-1) / 2


def euclid_dist(v):
    return (v[0] - v[1])**2


def per_class_recall(labels):
    recall_fxns = []
    for label_key in labels:
        label_idx = labels[label_key]
        fxn_name = label_key.replace('-', '_').replace(' ', '_')
        string_fxn = 'def ' + fxn_name + '_recall(y_true, y_pred):\n'
        string_fxn += '\ttrue_positives = K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0)\n'
        string_fxn += '\tpossible_positives = K.sum(K.round(K.clip(y_true, 0, 1)), axis=0)\n'
        string_fxn += '\treturn true_positives['+str(label_idx)+'] / (possible_positives['+str(label_idx)+'] + K.epsilon())\n'
          
        exec(string_fxn)
        recall_fxn = eval(fxn_name + '_recall')
        recall_fxns.append(recall_fxn)
    
    return recall_fxns


def per_class_precision(labels):
    precision_fxns = []
    
    for label_key in labels:
        label_idx = labels[label_key]
        fxn_name = label_key.replace('-', '_').replace(' ', '_')
        string_fxn = 'def ' + fxn_name + '_precision(y_true, y_pred):\n'
        string_fxn += '\ttrue_positives = K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0)\n'
        string_fxn += '\tpredicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)), axis=0)\n'
        string_fxn += '\treturn true_positives['+str(label_idx)+'] / (predicted_positives['+str(label_idx)+'] + K.epsilon())\n'
        
        exec(string_fxn)
        precision_fxn = eval(fxn_name + '_precision')
        precision_fxns.append(precision_fxn)
    
    return precision_fxns


def per_class_recall_3d(labels):
    recall_fxns = []
    
    for label_key in labels:
        label_idx = labels[label_key]
        fxn_prefix = label_key.replace('-', '_').replace(' ', '_')
        string_fxn = 'def ' + fxn_prefix + '_recall(y_true, y_pred):\n'
        string_fxn += '\ttrue_positives = K.sum(K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0), axis=0)\n'
        string_fxn += '\tpossible_positives = K.sum(K.sum(K.round(K.clip(y_true, 0, 1)), axis=0), axis=0)\n'
        string_fxn += '\treturn true_positives['+str(label_idx)+'] / (possible_positives['+str(label_idx)+'] + K.epsilon())\n'
        
        exec(string_fxn)
        recall_fxn = eval(fxn_prefix + '_recall')
        recall_fxns.append(recall_fxn)
    
    return recall_fxns


def per_class_precision_3d(labels):
    precision_fxns = []
    
    for label_key in labels:
        label_idx = labels[label_key]
        fxn_prefix = label_key.replace('-', '_').replace(' ', '_')
        string_fxn = 'def ' + fxn_prefix + '_precision(y_true, y_pred):\n'
        string_fxn += '\ttrue_positives = K.sum(K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0), axis=0)\n'
        string_fxn += '\tpredicted_positives = K.sum(K.sum(K.round(K.clip(y_pred, 0, 1)), axis=0), axis=0)\n'
        string_fxn += '\treturn true_positives['+str(label_idx)+'] / (predicted_positives['+str(label_idx)+'] + K.epsilon())\n'
        
        exec(string_fxn)
        precision_fxn = eval(fxn_prefix + '_precision')
        precision_fxns.append(precision_fxn)
    
    return precision_fxns


def per_class_recall_4d(labels):
    recall_fxns = []
    
    for label_key in labels:
        label_idx = labels[label_key]
        fxn_prefix = label_key.replace('-', '_').replace(' ', '_')
        string_fxn = 'def ' + fxn_prefix + '_recall(y_true, y_pred):\n'
        string_fxn += '\ttrue_positives = K.sum(K.sum(K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0), axis=0), axis=0)\n'
        string_fxn += '\tpossible_positives = K.sum(K.sum(K.sum(K.round(K.clip(y_true, 0, 1)), axis=0), axis=0), axis=0)\n'
        string_fxn += '\treturn true_positives['+str(label_idx)+'] / (possible_positives['+str(label_idx)+'] + K.epsilon())\n'
        
        exec(string_fxn)
        recall_fxn = eval(fxn_prefix + '_recall')
        recall_fxns.append(recall_fxn)
    
    return recall_fxns


def per_class_precision_4d(labels):
    precision_fxns = []
    
    for label_key in labels:
        label_idx = labels[label_key]
        fxn_prefix = label_key.replace('-', '_').replace(' ', '_')
        string_fxn = 'def ' + fxn_prefix + '_precision(y_true, y_pred):\n'
        string_fxn += '\ttrue_positives = K.sum(K.sum(K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0), axis=0), axis=0)\n'
        string_fxn += '\tpredicted_positives = K.sum(K.sum(K.sum(K.round(K.clip(y_pred, 0, 1)), axis=0), axis=0), axis=0)\n'
        string_fxn += '\treturn true_positives['+str(label_idx)+'] / (predicted_positives['+str(label_idx)+'] + K.epsilon())\n'
        
        exec(string_fxn)
        precision_fxn = eval(fxn_prefix + '_precision')
        precision_fxns.append(precision_fxn)
    
    return precision_fxns


def per_class_recall_5d(labels):
    recall_fxns = []
    
    for label_key in labels:
        label_idx = labels[label_key]
        fxn_prefix = label_key.replace('-', '_').replace(' ', '_')
        string_fxn = 'def ' + fxn_prefix + '_recall(y_true, y_pred):\n'
        string_fxn += '\ttrue_positives = K.sum(K.sum(K.sum(K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0), axis=0), axis=0), axis=0)\n'
        string_fxn += '\tpossible_positives = K.sum(K.sum(K.sum(K.sum(K.round(K.clip(y_true, 0, 1)), axis=0), axis=0), axis=0), axis=0)\n'
        string_fxn += '\treturn true_positives['+str(label_idx)+'] / (possible_positives['+str(label_idx)+'] + K.epsilon())\n'
        
        exec(string_fxn)
        recall_fxn = eval(fxn_prefix + '_recall')
        recall_fxns.append(recall_fxn)
    
    return recall_fxns


def per_class_precision_5d(labels):
    precision_fxns = []
    
    for label_key in labels:
        label_idx = labels[label_key]
        fxn_prefix = label_key.replace('-', '_').replace(' ', '_')
        string_fxn = 'def ' + fxn_prefix + '_precision(y_true, y_pred):\n'
        string_fxn += '\ttrue_positives = K.sum(K.sum(K.sum(K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0), axis=0), axis=0), axis=0)\n'
        string_fxn += '\tpredicted_positives = K.sum(K.sum(K.sum(K.sum(K.round(K.clip(y_pred, 0, 1)), axis=0), axis=0), axis=0), axis=0)\n'
        string_fxn += '\treturn true_positives['+str(label_idx)+'] / (predicted_positives['+str(label_idx)+'] + K.epsilon())\n'
        
        exec(string_fxn)
        precision_fxn = eval(fxn_prefix + '_precision')
        precision_fxns.append(precision_fxn)
    
    return precision_fxns


def get_metric_dict(output_tensor_maps):
    metrics = {}
    losses = []
    loss_weights = []
    for tm in output_tensor_maps:
        loss_weights.append(tm.loss_weight)
        for m in tm.metrics:
            if isinstance(m, str):
                metrics[m] = m
            else:
                metrics[m.__name__] = m
        if tm.loss == 'categorical_crossentropy':
            losses.append(categorical_crossentropy)
        elif tm.loss == 'binary_crossentropy':
            losses.append(binary_crossentropy)
        elif tm.loss == 'mean_absolute_error' or tm.loss == 'mae':
            losses.append(mean_absolute_error)
        elif tm.loss == 'mean_squared_error' or tm.loss == 'mse':
            losses.append(mean_squared_error)
        elif tm.loss == 'cosine_proximity':
            losses.append(cosine_proximity)
        elif tm.loss == 'logcosh':
            losses.append(logcosh)                        
        else:
            metrics[tm.loss.__name__] = tm.loss
            losses.append(tm.loss)

    def loss_fxn(y_true, y_pred):
        my_loss = 0
        for loss_fxn,loss_weight in zip(losses, loss_weights):
            my_loss += loss_weight*loss_fxn(y_true, y_pred)
        return my_loss
    metrics['loss'] = loss_fxn

    return metrics


def get_roc_aucs(predictions, truth, labels):
    """Compute ROC AUC for each label of each given model"""
    aucs = dict()

    for model in predictions.keys():
        roc_auc = dict()
        for label_name in labels.keys():
            label_encoding = labels[label_name]
            fpr, tpr,  _ = roc_curve(truth[:, label_encoding], predictions[model][:, label_encoding])
            roc_auc[label_name] = auc(fpr, tpr)
        aucs[model] = roc_auc

    return aucs


def get_precision_recall_aucs(predictions, truth, labels):
    """Compute Precision-Recall AUC for each label of each given model"""
    aucs = dict()

    for model in predictions.keys():
        average_precision = dict()
        for label_name in labels.keys():
            label_encoding = labels[label_name]
            average_precision[label_name] = average_precision_score(truth[:, label_encoding], predictions[model][:, label_encoding])
        aucs[model] = average_precision

    return aucs


def log_aucs(**aucs):
    """Log and tabulate AUCs given as nested dictionaries in the format '{model: {label: auc}}'"""
    def dashes(n): return '-' * n

    header = "{:<35} {:<20} {:<15}"
    row = "{:<35} {:<20} {:<15.10f}"
    width = 85
    logging.info(dashes(width))
    for auc_name, auc_value in aucs.items():
        logging.info(header.format('Model', 'Label', auc_name+' AUC'))
        for model, model_value in auc_value.items():
            for label, auc in model_value.items():
                logging.info(row.format(model, label, auc))
        logging.info(dashes(width))


def get_pearson_coefficients(predictions, truth):
    coefs = dict()
    for model in predictions.keys():
        # corrcoef() returns full covariance matrix
        pearson = np.corrcoef(predictions[model].flatten(), truth.flatten())[1, 0]
        coefs[model] = pearson

    return coefs


def log_pearson_coefficients(coefs, label):
    def dashes(n): return '-' * n

    header = "{:<30} {:<25} {:<15} {:<15}"
    row = "{:<30} {:<25} {:<15.10f} {:<15.10f}"
    width = 85
    logging.info(dashes(width))
    logging.info(header.format('Model', 'Label', 'Pearson R', 'Pearson R^2'))
    for model, coef in coefs.items():
        logging.info(row.format(model, label, coef, coef*coef))
    logging.info(dashes(width))
