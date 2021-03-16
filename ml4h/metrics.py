# metrics.py
import logging
import numpy as np
import tensorflow as tf
import tensorflow.keras.backend as K

from sklearn.metrics import roc_curve, auc, average_precision_score


from tensorflow.keras.losses import binary_crossentropy, categorical_crossentropy, logcosh, cosine_similarity, mean_squared_error, mean_absolute_error, mean_absolute_percentage_error

STRING_METRICS = [
    'categorical_crossentropy','binary_crossentropy','mean_absolute_error','mae',
    'mean_squared_error', 'mse', 'cosine_similarity', 'logcosh',
]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~ Metrics ~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def weighted_crossentropy(weights, name='anonymous'):
    """A weighted version of tensorflow.keras.objectives.categorical_crossentropy

    Arguments:
        weights = np.array([0.5,2,10]) # Class one at 0.5, class 2 twice the normal weights, class 3 10x.
        name: string identifying the loss to differentiate when models have multiple losses

    Returns:
        keras loss function named name+'_weighted_loss'

    """
    string_globe = 'global ' + name + '_weights\n'
    string_globe += 'global ' + name + '_kweights\n'
    string_globe += name + '_weights = np.array(weights)\n'
    string_globe += name + '_kweights = K.variable('+name+'_weights)\n'
    exec(string_globe, globals(), locals())
    fxn_postfix = '_weighted_loss'
    string_fxn = 'def ' + name + fxn_postfix + '(y_true, y_pred):\n'
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


def y_true_times_mse(y_true, y_pred):
    return K.maximum(y_true, 1.0)*mean_squared_error(y_true, y_pred)


def y_true_squared_times_mse(y_true, y_pred):
    return K.maximum(1.0+y_true, 1.0)*K.maximum(1.0+y_true, 1.0)*mean_squared_error(y_true, y_pred)


def y_true_cubed_times_mse(y_true, y_pred):
    return K.maximum(y_true, 1.0)*K.maximum(y_true, 1.0)*K.maximum(y_true, 1.0)*mean_squared_error(y_true, y_pred)


def y_true_squared_times_logcosh(y_true, y_pred):
    return K.maximum(1.0+y_true, 1.0)*K.maximum(1.0+y_true, 1.0)*logcosh(y_true, y_pred)


def asymmetric_outlier_mse(y_true, y_pred):
    """Loss function which asymmetrically penalizes over estimations of large values."""
    top_over = 40.0 * K.maximum(y_true - 2.0, 0.0) * K.maximum(y_true - y_pred, 0.0) * mean_squared_error(y_true, y_pred)
    top_over += 20.0 * K.maximum(y_true - 1.0, 0.0) * K.maximum(y_true - y_pred, 0.0) * mean_squared_error(y_true, y_pred)
    top_under = 5.0 * K.maximum(y_true - 1.0, 0.0) * K.maximum(y_pred - y_true, 0.0) * mean_squared_error(y_true, y_pred)
    return top_over + top_under + logcosh(y_true, y_pred)


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


def pearson(y_true, y_pred):
    # normalizing stage - setting a 0 mean.
    y_true -= K.mean(y_true, axis=-1)
    y_pred -= K.mean(y_pred, axis=-1)
    # normalizing stage - setting a 1 variance
    y_true = K.l2_normalize(y_true, axis=-1)
    y_pred = K.l2_normalize(y_pred, axis=-1)
    # final result
    pearson_correlation = K.sum(y_true * y_pred, axis=-1)
    return pearson_correlation


def _make_riskset(follow_up_times):
    # sort in descending order
    follow_up_times_np = tf.make_ndarray(tf.make_tensor_proto(follow_up_times))
    o = np.argsort(-follow_up_times_np)
    n_samples = follow_up_times_np.shape[0]
    risk_set = np.zeros((n_samples, n_samples))

    for i_start, i_sort in enumerate(o):
        time_i_start = follow_up_times_np[i_sort]
        k = i_start
        while k < n_samples and time_i_start <= follow_up_times_np[o[k]]:
            k += 1
        risk_set[i_sort, o[:k]] = True
    risk_set_tf = tf.convert_to_tensor(risk_set)
    return risk_set_tf


def _softmax_masked(risk_scores, mask, axis=0, keepdims=None):
    """Compute logsumexp across `axis` for entries where `mask` is true."""
    mask_f = K.cast(mask, risk_scores.dtype)
    risk_scores_masked = risk_scores * mask_f
    # for numerical stability, subtract the maximum value before taking the exponential
    amax = K.max(risk_scores_masked, axis=axis, keepdims=True)
    risk_scores_shift = risk_scores_masked - amax

    exp_masked = K.exp(risk_scores_shift) * mask_f
    exp_sum = K.sum(exp_masked, axis=axis, keepdims=True)
    output = amax + K.log(exp_sum)
    if not keepdims:
        output = K.squeeze(output, axis=axis)
    return output


@tf.function
def cox_hazard_loss(y_true, y_pred):
    # move batch dimension to the end so predictions get broadcast row-wise when multiplying by riskset
    pred_t = K.transpose(y_pred[:, 0])
    events = y_true[:, 0]
    follow_up_times = y_true[:, 1]
    # compute log of sum over risk set for each row
    rr = _softmax_masked(pred_t, _make_riskset(follow_up_times), axis=1, keepdims=True)

    losses = events * (rr - y_pred[:, 0])
    loss = K.mean(losses)
    return loss


def survival_likelihood_loss(n_intervals):
    """Create custom Keras loss function for neural network survival model.

    This function is tightly coupled with the function _survival_tensor defined in tensor_from_file.py which builds the y_true tensor.

    Arguments
        n_intervals: the number of survival time intervals
    Returns
        Custom loss function that can be used with Keras
    """

    def loss(y_true, y_pred):
        """
        To play nicely with the Keras framework y_pred is the same shape as y_true.
        However we only consider the first half (n_intervals) of y_pred.
        Arguments
            y_true: Tensor.
              First half of the values are 1 if individual survived that interval, 0 if not.
              Second half of the values are for individuals who failed, and are 1 for time interval during which failure occurred, 0 for other intervals.
              For example given n_intervals = 3 a sample with prevalent disease will have y_true [0, 0, 0, 1, 0, 0]
              a sample with incident disease occurring in the last time bin will have y_true [1, 1, 0, 0, 0, 1]
              a sample who is lost to follow up (censored) in middle time bin will have y_true [1, 0, 0, 0, 0, 0]
            y_pred: Tensor, predicted survival probability (1-hazard probability) for each time interval.
        Returns
            Vector of losses for this minibatch.
        """
        failure_likelihood = 1. - (y_true[:, n_intervals:] * y_pred[:, 0:n_intervals])  # Loss only for individuals who failed
        survival_likelihood = y_true[:, 0:n_intervals] * y_pred[:, 0:n_intervals]  # Loss for intervals that were survived
        survival_likelihood += 1. - y_true[:, 0:n_intervals]  # No survival loss if interval was censored or failed
        return K.sum(-K.log(K.clip(K.concatenate((survival_likelihood, failure_likelihood)), K.epsilon(), None)), axis=-1)  # return -log likelihood

    return loss


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
        elif tm.loss == 'cosine_similarity':
            losses.append(cosine_similarity)
        elif tm.loss == 'logcosh':
            losses.append(logcosh)
        elif tm.loss == 'mape':
            losses.append(mean_absolute_percentage_error)
        else:
            metrics[tm.loss.__name__] = tm.loss
            losses.append(tm.loss)

    def loss(y_true, y_pred):
        my_loss = 0
        for my_loss_fxn, loss_weight in zip(losses, loss_weights):
            my_loss += loss_weight * my_loss_fxn(y_true, y_pred)
        return my_loss
    metrics['loss'] = loss

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


def coefficient_of_determination(truth, predictions, eps=1e-6):
    true_mean = np.mean(truth)
    total_sum_of_squares = np.sum((truth - true_mean) * (truth - true_mean))
    residual_sum_of_squares = np.sum((predictions - truth) * (predictions - truth))
    r_squared = 1 - (residual_sum_of_squares / (total_sum_of_squares + eps))
    return r_squared


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


def _unpack_truth_into_events(truth, intervals):
    event_time = np.argmin(np.diff(truth[:, :intervals]), axis=-1)
    event_time[truth[:, intervals-1] == 1] = intervals-1  # If the sample is never censored set event time to max time
    event_indicator = np.sum(truth[:, intervals:], axis=-1).astype(np.bool)
    return event_indicator, event_time


def _get_comparable(event_indicator, event_time, order):
    n_samples = len(event_time)
    tied_time = 0
    comparable = {}
    i = 0
    while i < n_samples - 1:
        time_i = event_time[order[i]]
        start = i + 1
        end = start
        while end < n_samples and event_time[order[end]] == time_i:
            end += 1

        # check for tied event times
        event_at_same_time = event_indicator[order[i:end]]
        censored_at_same_time = ~event_at_same_time
        for j in range(i, end):
            if event_indicator[order[j]]:
                mask = np.zeros(n_samples, dtype=bool)
                mask[end:] = True
                # an event is comparable to censored samples at same time point
                mask[i:end] = censored_at_same_time
                comparable[j] = mask
                tied_time += censored_at_same_time.sum()
        i = end

    return comparable, tied_time


def concordance_index(prediction, truth, tied_tol=1e-8):
    intervals = truth.shape[-1] // 2
    event_indicator, event_time = _unpack_truth_into_events(truth, intervals)
    estimate = np.cumprod(prediction[:, :intervals], axis=-1)[:, -1]
    order = np.argsort(event_time)
    comparable, tied_time = _get_comparable(event_indicator, event_time, order)

    concordant = 0
    discordant = 0
    tied_risk = 0
    numerator = 0.0
    denominator = 0.0
    for ind, mask in comparable.items():
        est_i = estimate[order[ind]]
        event_i = event_indicator[order[ind]]

        est = estimate[order[mask]]

        assert event_i, 'got censored sample at index %d, but expected uncensored' % order[ind]

        ties = np.absolute(est - est_i) <= tied_tol
        n_ties = ties.sum()

        # an event should have a higher score
        con = est > est_i
        n_con = con[~ties].sum()

        numerator += n_con + 0.5 * n_ties
        denominator += mask.sum()

        tied_risk += n_ties
        concordant += n_con
        discordant += est.size - n_con - n_ties

    cindex = numerator / denominator
    return cindex, concordant, discordant, tied_risk, tied_time
