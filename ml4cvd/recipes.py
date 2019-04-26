# recipes.py

# Imports
import os
import csv
import logging
import numpy as np
from functools import reduce
from timeit import default_timer as timer
from collections import Counter, defaultdict

from arguments import parse_args
from defines import TENSOR_EXT
from tensor_writer_ukbb import write_tensors
from tensor_map_maker import write_tensor_maps
from plots import evaluate_predictions, plot_scatters, plot_rocs, plot_precision_recalls
from explorations import sample_from_char_model, mri_dates, ecg_dates, predictions_to_pngs, plot_while_learning
from metrics import get_roc_aucs, get_precision_recall_aucs, get_pearson_coefficients, log_aucs, log_pearson_coefficients
from tensor_generators import TensorGenerator, test_train_valid_tensor_generators, big_batch_from_minibatch_generator, get_test_train_valid_paths
from models import make_multimodal_to_multilabel_model, train_model_from_generators, get_model_inputs_outputs, make_shallow_model, make_character_model_plus


def run(args):
    try:
        # Keep track of elapsed execution time
        start_time = timer()

        if 'tensorize' == args.mode:
            write_tensors(args.id, args.db, args.xml_folder, args.zip_folder, args.phenos_folder, args.output_folder,
                          args.tensors, args.dicoms, args.volume_csv, args.lv_mass_csv, args.icd_csv, args.categorical_field_ids,
                          args.continuous_field_ids, args.mri_field_ids, args.xml_field_ids, args.x, args.y, args.z,
                          args.include_heart_zoom, args.zoom_x, args.zoom_y, args.zoom_width,  args.zoom_height,
                          args.write_pngs, args.min_sample_id, args.max_sample_id, args.min_values)
        elif 'train' == args.mode:
            train_multimodal_multitask(args)
        elif 'test' == args.mode:
            test_multimodal_multitask(args)
        elif 'compare' == args.mode:
            compare_multimodal_multitask_models(args)
        elif 'infer' == args.mode:
            infer_multimodal_multitask(args)
        elif 'segmentation_to_pngs' == args.mode:
            segmentation_to_pngs(args)
        elif 'plot_while_training' == args.mode:
            plot_while_training(args)
        elif 'plot_mri_dates' == args.mode:
            mri_dates(args.tensors, args.output_folder, args.id)
        elif 'plot_ecg_dates' == args.mode:
            ecg_dates(args.tensors, args.output_folder, args.id)
        elif 'train_shallow' == args.mode:
            train_shallow_model(args)
        elif 'train_char' == args.mode:
            train_char_model(args)
        elif 'write_tensor_maps' == args.mode:
            write_tensor_maps(args)
        else:
            raise ValueError('Unknown mode:', args.mode)

    except Exception as e:
        logging.exception(e)

    end_time = timer()
    elapsed_time = end_time - start_time
    logging.info("Executed the '{}' operation in {:.2f} seconds".format(args.mode, elapsed_time))


def train_multimodal_multitask(args):
    generate_train, generate_valid, generate_test = test_train_valid_tensor_generators(args.tensor_maps_in, args.tensor_maps_out, args.tensors,
                                                                                       args.batch_size, args.valid_ratio, args.test_ratio,
                                                                                       args.icd_csv, args.balance_by_icds)

    model = make_multimodal_to_multilabel_model(args.model_file, args.model_layers, args.model_freeze, args.tensor_maps_in, args.tensor_maps_out,
                                                args.activation, args.dense_layers, args.dropout, args.mlp_concat, args.conv_layers, args.max_pools,
                                                args.res_layers, args.dense_blocks, args.block_size, args.conv_bn, args.conv_x, args.conv_y,
                                                args.conv_z, args.conv_dropout, args.conv_width, args.u_connect, args.pool_x, args.pool_y,
                                                args.pool_z, args.padding, args.learning_rate)

    model = train_model_from_generators(model, generate_train, generate_valid, args.training_steps, args.validation_steps, args.batch_size,
                                        args.epochs, args.patience, args.output_folder, args.id, args.inspect_model, args.inspect_show_labels)

    test_data, test_labels, test_paths = big_batch_from_minibatch_generator(args.tensor_maps_in, args.tensor_maps_out, generate_test, args.test_steps)
    return _predict_and_evaluate(model, test_data, test_labels, args.tensor_maps_out, args.batch_size, args.output_folder, args.id, test_paths)


def test_multimodal_multitask(args):
    _, _, generate_test = test_train_valid_tensor_generators(args.tensor_maps_in, args.tensor_maps_out, args.tensors, args.batch_size,
                                                             args.valid_ratio, args.test_ratio, args.icd_csv, args.balance_by_icds)

    model = make_multimodal_to_multilabel_model(args.model_file, args.model_layers, args.model_freeze, args.tensor_maps_in, args.tensor_maps_out,
                                                args.activation, args.dense_layers, args.dropout, args.mlp_concat, args.conv_layers, args.max_pools,
                                                args.res_layers, args.dense_blocks, args.block_size, args.conv_bn, args.conv_x, args.conv_y,
                                                args.conv_z, args.conv_dropout, args.conv_width, args.u_connect, args.pool_x, args.pool_y,
                                                args.pool_z, args.padding, args.learning_rate)

    data, labels, paths = big_batch_from_minibatch_generator(args.tensor_maps_in, args.tensor_maps_out, generate_test, args.test_steps)
    return _predict_and_evaluate(model, data, labels, args.tensor_maps_out, args.batch_size, args.output_folder, args.id, paths)
  

def compare_multimodal_multitask_models(args):
    input_prefix = "input"
    output_prefix = "output"

    tensor_paths = _get_tensor_files(args.tensors)
    generator = TensorGenerator(args.batch_size, args.tensor_maps_in, args.tensor_maps_out, tensor_paths, keep_paths=True)

    input_data, output_data, paths = big_batch_from_minibatch_generator(args.tensor_maps_in, args.tensor_maps_out, generator, args.test_steps)
    models_inputs_outputs = get_model_inputs_outputs(args.model_files, args.tensor_maps_in, args.tensor_maps_out)
    common_outputs = _get_common_outputs(models_inputs_outputs, output_prefix)
    predictions = _get_predictions(args, models_inputs_outputs, input_data, common_outputs, input_prefix, output_prefix)
    _calculate_and_plot_prediction_stats(args, predictions, output_data)


def infer_multimodal_multitask(args):
    stats = Counter()
    tensor_paths_inferred = {}
    tensor_paths = [ args.tensors + tp for tp in os.listdir(args.tensors) if os.path.splitext(tp)[-1].lower()==TENSOR_EXT ]
    # hard code batch size to 1 so we can iterate over filenames and generated tensors together in the tensor_paths for loop
    generate_test = TensorGenerator(1, args.tensor_maps_in, args.tensor_maps_out, tensor_paths, keep_paths=True)
    model = make_multimodal_to_multilabel_model(args.model_file, args.model_layers, args.model_freeze, args.tensor_maps_in, args.tensor_maps_out,
                                                args.activation, args.dense_layers, args.dropout, args.mlp_concat, args.conv_layers, args.max_pools,
                                                args.res_layers, args.dense_blocks, args.block_size, args.conv_bn, args.conv_x, args.conv_y,
                                                args.conv_z, args.conv_dropout, args.conv_width, args.u_connect, args.pool_x, args.pool_y,
                                                args.pool_z, args.padding, args.learning_rate)
    
    with open(os.path.join(args.output_folder, args.id, 'inference_' + args.id + '.tsv' ), mode='w') as inference_file:
        inference_writer = csv.writer(inference_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        inference_writer.writerow(['sample_id'] + [ot for ot,otm in zip(args.output_tensors, args.tensor_maps_out) if len(otm.shape)==1])

        while True:
            tensor_batch = next(generate_test)
            tensor_path = tensor_batch[2][0]
            if tensor_path in tensor_paths_inferred:
                print('Done inferring tensors. Looped over at tensor:', tensor_path)
                break

            prediction = model.predict(tensor_batch[0])
            if len(args.tensor_maps_out) == 1:
                prediction = [prediction]

            csv_row = [os.path.basename(tensor_path).replace(TENSOR_EXT, '')] # extract sample id
            for y,tm in zip(prediction, args.tensor_maps_out):
                if len(tm.shape) == 1:
                    csv_row.append(str(tm.rescale(y)[0][0]))  # first index into batch then index into the 1x1 structure
            inference_writer.writerow(csv_row)

            tensor_paths_inferred[tensor_path] = True
            stats['count'] += 1
            if stats['count']% 500 == 0:
                print('Wrote:', stats['count'], 'rows of inference.  Last tensor:', tensor_path)


def train_shallow_model(args):
    generate_train, generate_valid, generate_test = test_train_valid_tensor_generators(args.tensor_maps_in, args.tensor_maps_out, args.tensors,
                                                                                       args.batch_size, args.valid_ratio, args.test_ratio,
                                                                                       args.icd_csv, args.balance_by_icds)

    model = make_shallow_model(args.tensor_maps_in, args.tensor_maps_out, args.learning_rate, args.model_file, args.model_layers)
    model = train_model_from_generators(model, generate_train, generate_valid, args.training_steps, args.validation_steps, args.batch_size,
                                        args.epochs, args.patience, args.output_folder, args.id, args.inspect_model, args.inspect_show_labels)

    test_data, test_labels, test_paths = big_batch_from_minibatch_generator(args.tensor_maps_in, args.tensor_maps_out, generate_test, args.test_steps)
    return _predict_and_evaluate(model, test_data, test_labels, args.tensor_maps_out, args.batch_size, args.output_folder, args.id, test_paths)


def train_char_model(args):
    base_model = make_multimodal_to_multilabel_model(args.model_file, args.model_layers, args.model_freeze, args.tensor_maps_in, args.tensor_maps_out,
                                                     args.activation, args.dense_layers, args.dropout, args.mlp_concat, args.conv_layers,
                                                     args.max_pools, args.res_layers, args.dense_blocks, args.block_size, args.conv_bn, args.conv_x,
                                                     args.conv_y, args.conv_z, args.conv_dropout, args.conv_width, args.u_connect, args.pool_x,
                                                     args.pool_y, args.pool_z, args.padding, args.learning_rate)

    model, char_model = make_character_model_plus(args.tensor_maps_in, args.tensor_maps_out, args.learning_rate, base_model, args.model_layers)
    generate_train, generate_valid, generate_test = test_train_valid_tensor_generators(args.tensor_maps_in, args.tensor_maps_out, args.tensors,
                                                                                       args.batch_size, args.valid_ratio, args.test_ratio,
                                                                                       args.icd_csv, args.balance_by_icds)

    model = train_model_from_generators(model, generate_train, generate_valid, args.training_steps, args.validation_steps, args.batch_size,
                                        args.epochs, args.patience, args.output_folder, args.id, args.inspect_model, args.inspect_show_labels)
    test_batch, _, test_paths = next(generate_test)
    sample_from_char_model(char_model, test_batch, test_paths)
    data, labels, paths = big_batch_from_minibatch_generator(args.tensor_maps_in, args.tensor_maps_out, generate_test, args.test_steps)
    return _predict_and_evaluate(model, data, labels, args.tensor_maps_out, args.batch_size, args.output_folder, args.id, paths)


def segmentation_to_pngs(args):
    _, _, generate_test = test_train_valid_tensor_generators(args.tensor_maps_in, args.tensor_maps_out, args.tensors, args.batch_size*args.test_steps,
                                                             args.valid_ratio, args.test_ratio, args.icd_csv, args.balance_by_icds, True)

    model = make_multimodal_to_multilabel_model(args.model_file, args.model_layers, args.model_freeze, args.tensor_maps_in, args.tensor_maps_out,
                                                args.activation, args.dense_layers, args.dropout, args.mlp_concat, args.conv_layers, args.max_pools,
                                                args.res_layers, args.dense_blocks, args.block_size, args.conv_bn, args.conv_x, args.conv_y,
                                                args.conv_z, args.conv_dropout, args.conv_width, args.u_connect, args.pool_x, args.pool_y,
                                                args.pool_z, args.padding, args.learning_rate)

    data, labels, paths = big_batch_from_minibatch_generator(args.tensor_maps_in, args.tensor_maps_out, generate_test, args.test_steps)
    predictions = model.predict(data, batch_size=args.batch_size)
    if len(args.tensor_maps_out) == 1:
        predictions = [predictions]
    folder = os.path.join(args.output_folder, args.id) + '/'
    predictions_to_pngs(predictions, args.tensor_maps_in, args.tensor_maps_out, data, labels, paths, folder)


def plot_while_training(args):
    generate_train, _, generate_test = test_train_valid_tensor_generators(args.tensor_maps_in, args.tensor_maps_out, args.tensors, args.batch_size,
                                                                          args.valid_ratio, args.test_ratio, args.icd_csv, args.balance_by_icds)

    test_data, test_labels, test_paths = big_batch_from_minibatch_generator(args.tensor_maps_in, args.tensor_maps_out, generate_test, args.test_steps)
    model = make_multimodal_to_multilabel_model(args.model_file, args.model_layers, args.model_freeze, args.tensor_maps_in, args.tensor_maps_out,
                                                args.activation, args.dense_layers, args.dropout, args.mlp_concat, args.conv_layers, args.max_pools,
                                                args.res_layers, args.dense_blocks, args.block_size, args.conv_bn, args.conv_x, args.conv_y,
                                                args.conv_z, args.conv_dropout, args.conv_width, args.u_connect, args.pool_x, args.pool_y,
                                                args.pool_z, args.padding, args.learning_rate)

    plot_folder = os.path.join(args.output_folder, args.id, 'training_frames/')
    plot_while_learning(model, args.tensor_maps_in, args.tensor_maps_out, generate_train, test_data, test_labels, test_paths, args.epochs,
                        args.batch_size, args.training_steps, plot_folder, args.id, args.write_pngs)


def _predict_and_evaluate(model, test_data, test_labels, tensor_maps_out, batch_size, output_folder, run_id, test_paths=None):
    performance_metrics = {}
    plot_path = os.path.join(output_folder, run_id)
    y_pred = model.predict(test_data, batch_size=batch_size)
    for y, tm in zip(y_pred, tensor_maps_out):
        if len(tensor_maps_out) == 1:
            y = y_pred
        performance_metrics.update(evaluate_predictions(tm, y, test_labels, test_data, tm.name, plot_path, test_paths))
    return performance_metrics


def _get_common_outputs(models_inputs_outputs, output_prefix):
    """Returns a set of outputs common to all the models so we can compare the models according to those outputs only"""
    all_outputs = []
    for (_, ios) in models_inputs_outputs.items():
        outputs = {k: v for (k, v) in ios.items() if k == output_prefix}
        for (_, output) in outputs.items():
            all_outputs.append(set(output))
    return reduce(set.intersection, all_outputs)


def _get_predictions(args, models_inputs_outputs, input_data, outputs, input_prefix, output_prefix):
    """Makes multi-modal predictions for a given number of models.

    Returns:
        dict: The nested dictionary of predicted values.

            {
                'tensor_map_1':
                    {
                        'model_1': [[value1, value2], [...]],
                        'model_2': [[value3, value4], [...]]
                    },
                'tensor_map_2':
                    {
                        'model_2': [[value5, value6], [...]],
                        'model_3': [[value7, value8], [...]]
                    }
            }
    """
    predictions = defaultdict(dict)
    for model_file in models_inputs_outputs.keys():
        args.model_file = model_file
        args.tensor_maps_in = models_inputs_outputs[model_file][input_prefix]
        args.tensor_maps_out = models_inputs_outputs[model_file][output_prefix]

        model = make_multimodal_to_multilabel_model(args.model_file, args.model_layers, args.model_freeze,
                                                    args.tensor_maps_in, args.tensor_maps_out, args.activation,
                                                    args.dense_layers, args.dropout, args.mlp_concat, args.conv_layers,
                                                    args.max_pools, args.res_layers, args.dense_blocks, args.block_size,
                                                    args.conv_bn, args.conv_x, args.conv_y, args.conv_z,
                                                    args.conv_dropout, args.conv_width, args.u_connect, args.pool_z,
                                                    args.padding, args.learning_rate)

        model_name = os.path.basename(model_file).replace(TENSOR_EXT, '')

        # We can feed 'model.predict()' the entire input data because it knows what subset to use
        y_pred = model.predict(input_data, batch_size=args.batch_size)

        for i, tm in enumerate(args.tensor_maps_out):
            if tm in outputs:
                if len(args.tensor_maps_out) == 1:
                    predictions[tm][model_name] = y_pred
                else:
                    predictions[tm][model_name] = y_pred[i]

    return predictions


def _calculate_and_plot_prediction_stats(args, predictions, outputs):
    for tm in args.tensor_maps_out:
        plot_title = tm.name+'_'+args.id
        plot_folder = os.path.join(args.output_folder, args.id)

        if tm.is_categorical_any_with_shape_len(1):
            msg = "For tm '{}' with channel map {}: sum truth = {}; sum pred = {}"
            for m in predictions[tm]:
                logging.info(msg.format(tm.name, tm.channel_map, np.sum(outputs[tm.output_name()], axis=0), np.sum(predictions[tm][m], axis=0)))
            plot_rocs(predictions[tm], outputs[tm.output_name()], tm.channel_map, plot_title, plot_folder)
        elif tm.is_categorical_any_with_shape_len(4):
            for p in predictions[tm]:
                y = predictions[tm][p]
                melt_shape = (y.shape[0]*y.shape[1]*y.shape[2]*y.shape[3], y.shape[4])
                predictions[tm][p] = y.reshape(melt_shape)

            y_truth = outputs[tm.output_name()].reshape(melt_shape)
            plot_rocs(predictions[tm], y_truth, tm.channel_map, plot_title, plot_folder)
            plot_precision_recalls(predictions[tm], y_truth, tm.channel_map, plot_title, plot_folder)

            roc_aucs = get_roc_aucs(predictions[tm], y_truth, tm.channel_map)
            precision_recall_aucs = get_precision_recall_aucs(predictions[tm], y_truth, tm.channel_map)

            aucs = {"ROC": roc_aucs, "Precision-Recall": precision_recall_aucs}
            log_aucs(**aucs)
        else:
            plot_scatters(predictions[tm], tm.rescale(outputs[tm.output_name()]), plot_title, plot_folder)

            coefs = get_pearson_coefficients(predictions[tm], tm.rescale(outputs[tm.output_name()]))
            log_pearson_coefficients(coefs, tm.name)


def _get_tensor_files(tensor_dir):
    return [tensor_dir + tp for tp in os.listdir(args.tensors) if os.path.splitext(tp)[-1].lower() == TENSOR_EXT]


if __name__=='__main__':
    args = parse_args()
    run(args)  # back to the top
