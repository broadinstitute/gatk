# recipes.py

# Imports
import os
import csv
from typing import Dict, List
from operator import itemgetter

import h5py
import copy
import logging
import numpy as np
import pandas as pd
from functools import reduce
from timeit import default_timer as timer
from collections import Counter, defaultdict

from ml4cvd.arguments import parse_args
from ml4cvd.TensorMap import Interpretation
from ml4cvd.defines import TENSOR_EXT, MODEL_EXT
from ml4cvd.tensor_map_maker import write_tensor_maps
from ml4cvd.explorations import sample_from_char_model, mri_dates, ecg_dates, predictions_to_pngs, sort_csv
from ml4cvd.explorations import tabulate_correlations_of_tensors, test_labels_to_label_map, infer_with_pixels
from ml4cvd.explorations import plot_heatmap_of_tensors, plot_while_learning, plot_histograms_of_tensors_in_pdf
from ml4cvd.tensor_generators import TensorGenerator, test_train_valid_tensor_generators, big_batch_from_minibatch_generator
from ml4cvd.models import make_character_model_plus, embed_model_predict, make_siamese_model, make_multimodal_multitask_model
from ml4cvd.plots import evaluate_predictions, plot_scatters, plot_rocs, plot_precision_recalls, plot_roc_per_class, plot_tsne
from ml4cvd.metrics import get_roc_aucs, get_precision_recall_aucs, get_pearson_coefficients, log_aucs, log_pearson_coefficients
from ml4cvd.plots import subplot_rocs, subplot_comparison_rocs, subplot_scatters, subplot_comparison_scatters, plot_saliency_maps, plot_partners_ecgs
from ml4cvd.tensor_writer_ukbb import write_tensors, append_fields_from_csv, append_gene_csv, write_tensors_from_dicom_pngs, write_tensors_from_ecg_pngs
from ml4cvd.models import train_model_from_generators, get_model_inputs_outputs, make_shallow_model, make_hidden_layer_model, saliency_map, make_variational_multimodal_multitask_model


def run(args):
    try:
        # Keep track of elapsed execution time
        start_time = timer()

        if 'tensorize' == args.mode:
            write_tensors(args.id, args.xml_folder, args.zip_folder, args.output_folder, args.tensors, args.dicoms, args.mri_field_ids, args.xml_field_ids,
                          args.zoom_x, args.zoom_y, args.zoom_width, args.zoom_height, args.write_pngs, args.min_sample_id,
                          args.max_sample_id, args.min_values)
        elif 'tensorize_pngs' == args.mode:
            write_tensors_from_dicom_pngs(args.tensors, args.dicoms, args.app_csv, args.dicom_series, args.min_sample_id, args.max_sample_id)
        elif 'tensorize_ecg_pngs' == args.mode:
            write_tensors_from_ecg_pngs(args.tensors, args.xml_folder, args.min_sample_id, args.max_sample_id)
        elif 'explore' == args.mode:
            explore(args)
        elif 'train' == args.mode:
            train_multimodal_multitask(args)
        elif 'test' == args.mode:
            test_multimodal_multitask(args)
        elif 'compare' == args.mode:
            compare_multimodal_multitask_models(args)
        elif 'infer' == args.mode:
            infer_multimodal_multitask(args)
        elif 'infer_hidden' == args.mode:
            infer_hidden_layer_multimodal_multitask(args)
        elif 'infer_pixels' == args.mode:
            infer_with_pixels(args)
        elif 'test_scalar' == args.mode:
            test_multimodal_scalar_tasks(args)
        elif 'compare_scalar' == args.mode:
            compare_multimodal_scalar_task_models(args)
        elif 'plot_predictions' == args.mode:
            plot_predictions(args)
        elif 'plot_while_training' == args.mode:
            plot_while_training(args)
        elif 'plot_saliency' == args.mode:
            saliency_maps(args)
        elif 'plot_mri_dates' == args.mode:
            mri_dates(args.tensors, args.output_folder, args.id)
        elif 'plot_ecg_dates' == args.mode:
            ecg_dates(args.tensors, args.output_folder, args.id)
        elif 'plot_histograms' == args.mode:
            plot_histograms_of_tensors_in_pdf(args.id, args.tensors, args.output_folder, args.max_samples)
        elif 'plot_heatmap' == args.mode:
            plot_heatmap_of_tensors(args.id, args.tensors, args.output_folder, args.min_samples, args.max_samples)
        elif 'plot_partners_ecgs' == args.mode:
            plot_partners_ecgs(args)
        elif 'tabulate_correlations' == args.mode:
            tabulate_correlations_of_tensors(args.id, args.tensors, args.output_folder, args.min_samples, args.max_samples)
        elif 'train_shallow' == args.mode:
            train_shallow_model(args)
        elif 'train_char' == args.mode:
            train_char_model(args)
        elif 'train_siamese' == args.mode:
            train_siamese_model(args)
        elif 'write_tensor_maps' == args.mode:
            write_tensor_maps(args)
        elif 'sort_csv' == args.mode:
            sort_csv(args.tensors, args.tensor_maps_in)
        elif 'append_continuous_csv' == args.mode:
            append_fields_from_csv(args.tensors, args.app_csv, 'continuous', ',')
        elif 'append_categorical_csv' == args.mode:
            append_fields_from_csv(args.tensors, args.app_csv, 'categorical', ',')
        elif 'append_continuous_tsv' == args.mode:
            append_fields_from_csv(args.tensors, args.app_csv, 'continuous', '\t')
        elif 'append_categorical_tsv' == args.mode:
            append_fields_from_csv(args.tensors, args.app_csv, 'categorical', '\t')
        elif 'append_gene_csv' == args.mode:
            append_gene_csv(args.tensors, args.app_csv, ',')
        else:
            raise ValueError('Unknown mode:', args.mode)

    except Exception as e:
        logging.exception(e)

    end_time = timer()
    elapsed_time = end_time - start_time
    logging.info("Executed the '{}' operation in {:.2f} seconds".format(args.mode, elapsed_time))


def _init_dict_of_tensors(tmaps: list) -> dict:
    # Iterate through tmaps and initialize dict in which to store
    # tensors, error types, and fpath
    dict_of_tensors = defaultdict(dict)
    for tm in tmaps:
        if tm.channel_map:
            for cm in tm.channel_map:
                dict_of_tensors[tm.name].update({(tm.name, cm): list()})
        else:
            dict_of_tensors[tm.name].update({f"{tm.name}": list()})
        dict_of_tensors[tm.name].update({f"error_type_{tm.name}": list()})
        dict_of_tensors[tm.name].update({"fpath": list()})
    return dict_of_tensors


def _tensors_to_df(args):
    generators = test_train_valid_tensor_generators(**args.__dict__)
    tmaps = [tm for tm in args.tensor_maps_in]
    list_of_tensor_dicts: List[Dict] = []
    dependents = {}
    for gen in generators:
        data_len = len(gen.path_iters[0].paths)
        for i, path in enumerate(gen.path_iters[0].paths):
            if (i+1) % 500 == 0:
                logging.info(f"{gen.name} - Parsing {i}/{data_len} ({i/data_len*100:.1f}%) done")
            tensor_dict = _init_dict_of_tensors(tmaps)
            try:
                with h5py.File(path, "r") as hd5:
                    # Iterate through each tmap
                    for tm in tmaps:
                        error_type = ""
                        try:
                            tensor = tm.tensor_from_file(tm, hd5, dependents)
                            tensor = tm.postprocess_tensor(tensor, augment=False, hd5=hd5)
                           
                            # Append tensor to dict
                            if tm.channel_map:
                                for cm in tm.channel_map:
                                    tensor_dict[tm.name][(tm.name, cm)] = tensor[tm.channel_map[cm]]
                            else:
                                # If tensor is a scalar, isolate the value in the array;
                                # otherwise, retain the value as array
                                if tm.shape[0] == 1:
                                    tensor = tensor.item()
                                tensor_dict[tm.name][tm.name] = tensor
                        except (IndexError, KeyError, ValueError, OSError, RuntimeError) as e:
                            # Could not obtain tensor, so append nans
                            if tm.channel_map:
                                for cm in tm.channel_map:
                                    tensor_dict[tm.name][(tm.name, cm)] = np.nan
                            else:
                                # TODO figure out the np.full stuff
                                tensor_dict[tm.name][tm.name] = np.full(tm.shape, np.nan)[0]

                            # Save  error type to more readable string
                            error_type = type(e).__name__

                        # Save error type, fpath, and generator name (set)
                        tensor_dict[tm.name][f"error_type_{tm.name}"] = error_type
                        tensor_dict[tm.name]["fpath"] = path
                        tensor_dict[tm.name]["generator"] = gen.name
            except OSError as e:
                logging.info(f"OSError {e}")

            # Append list of dicts with tensor_dict
            list_of_tensor_dicts.append(tensor_dict)

    # Now we have a list of dicts where each dict has {tmaps:values} and
    # each HD5 -> one dict in the list
    # Next, convert list of dicts -> dataframe
    df = pd.DataFrame()

    for tm in tmaps:
        # Isolate all {tmap:values} from the list of dicts for this tmap
        list_of_tmap_dicts = list(map(itemgetter(tm.name), list_of_tensor_dicts))

        # Convert this tmap-specific list of dicts into dict of lists
        dict_of_tmap_lists = {k: [d[k] for d in list_of_tmap_dicts]
                                 for k in list_of_tmap_dicts[0]}

        # Convert list of dicts into dataframe and concatenate to big df
        df = pd.concat([df, pd.DataFrame(dict_of_tmap_lists)], axis=1)

    # Remove duplicate columns: error_types, fpath
    df = df.loc[:, ~df.columns.duplicated()]

    # Remove "_worker" from "generator" values
    df["generator"].replace("_worker", "", regex=True, inplace=True)

    # Rearrange df columns so fpath and generator are at the end
    cols = [col for col in df if col not in ["fpath", "generator"]] \
           + ["fpath", "generator"]
    df = df[cols]

    # Cast dtype of some columns to string. Note this is necessary; although a
    # df (or pd.series) of floats will have the type "float", a df of strings
    # assumes a dtype of "object". Casting to dtype "string" will confer performnace
    # improvements in future versions of Pandas
    df["fpath"] = df["fpath"].astype("string")
    df["generator"] = df["generator"].astype("string")

    # Iterate through tensor (and channel) maps and cast Pandas dtype to string
    if Interpretation.LANGUAGE in [tm.interpretation for tm in tmaps]:
        for tm in [tm for tm in args.tensor_maps_in if tm.interpretation is Interpretation.LANGUAGE]: 
            if tm.channel_map:
                for cm in tm.channel_map:
                    key = (tm.name, cm)
                    df[key] = df[key].astype("string")
            else:
                key = tm.name
                df[key] = df[key].astype("string")
    logging.info(f"Extracted {len(tmaps)} tmaps from {df.shape[0]} hd5 files into DataFrame")
    return df


def explore(args):
    args.num_workers = 0
    tmaps = args.tensor_maps_in
    fpath_prefix = "summary_stats"

    if any([len(tm.shape) != 1 for tm in tmaps]):
        raise ValueError("Explore only works for 1D tensor maps, but len(tm.shape) returned a value other than 1.")

    # Iterate through tensors, get tmaps, and save to dataframe
    df = _tensors_to_df(args)
    
    # Save dataframe to CSV
    fpath = os.path.join(args.output_folder, args.id, "tensors_all_union.csv")
    df.to_csv(fpath, index=False)
    fpath = os.path.join(args.output_folder, args.id, "tensors_all_intersect.csv")
    df.dropna().to_csv(fpath, index=False)
    logging.info(f"Saved dataframe of tensors (union and intersect) to {fpath}")

    #fpath = os.path.join(args.output_folder, args.id, "tensors_all_union.csv")
    #logging.info(f"Loading {fpath} to Pandas DataFrame")
    #df = pd.read_csv(fpath, keep_default_na=False)
    #logging.info(f"Loaded {fpath} to Pandas DataFrame")

    # Check if any tmaps are categorical
    if Interpretation.CATEGORICAL in [tm.interpretation for tm in tmaps]:

        # Iterate through 1) df, 2) df without NaN-containing rows (intersect)
        for df_cur, df_str in zip([df, df.dropna()], ["union", "intersect"]):
            for tm in [tm for tm in tmaps if tm.interpretation is Interpretation.CATEGORICAL]:
                counts = []
                counts_missing = []
                if tm.channel_map:
                    for cm in tm.channel_map:
                        key = (tm.name, cm)
                        counts.append(df_cur[key].sum())
                        counts_missing.append(df_cur[key].isna().sum())
                else:
                    key = tm.name
                    counts.append(df_cur[key].sum())
                    counts_missing.append(df_cur[key].isna().sum())

                # Append list with missing counts
                counts.append(counts_missing[0])

                # Append list with total counts
                counts.append(sum(counts))

                # Create list of row names
                cm_names = [cm for cm in tm.channel_map] + [f"missing", f"total"]

                # Transform list into dataframe indexed by channel maps
                df_stats = pd.DataFrame(counts, index=cm_names, columns=["counts"])

                # Add new column: percent of all counts
                df_stats["fraction_of_total"] = df_stats["counts"] / df_stats.loc[f"total"]["counts"]

                # Save parent dataframe to CSV on disk
                fpath = os.path.join(args.output_folder, args.id, 
                            f"{fpath_prefix}_{Interpretation.CATEGORICAL}_{tm.name}_{df_str}.csv")
                df_stats.to_csv(fpath)
                logging.info(f"Saved summary stats of {Interpretation.CATEGORICAL} {tm.name} tmaps to {fpath}")

    # Check if any tmaps are continuous
    if Interpretation.CONTINUOUS in [tm.interpretation for tm in tmaps]:

        # Iterate through 1) df, 2) df without NaN-containing rows (intersect)
        for df_cur, df_str in zip([df, df.dropna()], ["union", "intersect"]):
            df_stats = pd.DataFrame()
            if df_cur.empty:
                logging.info(f"{df_str} of tensors results in empty dataframe. Skipping calculations of {Interpretation.CONTINUOUS} summary statistics")
            else:
                for tm in [tm for tm in tmaps if tm.interpretation is Interpretation.CONTINUOUS]:
                    if tm.channel_map:
                        for cm in tm.channel_map:
                            stats = dict()
                            key = (tm.name, cm)
                            stats["min"] = df_cur[key].min()
                            stats["max"] = df_cur[key].max()
                            stats["mean"] = df_cur[key].mean()
                            stats["median"] = df_cur[key].median()
                            stats["mode"] = df_cur[key].mode()[0]
                            stats["variance"] = df_cur[key].var()
                            stats["count"] = df_cur[key].count()
                            stats["missing"] = df_cur[key].isna().sum()
                            stats["total"] = len(df_cur[key])
                            stats["missing_fraction"] = stats["missing"] / stats["total"]
                            df_stats = pd.concat([df_stats, pd.DataFrame([stats], index=[cm])])
                    else:
                        stats = dict()
                        key = tm.name
                        stats["min"] = df_cur[key].min()
                        stats["max"] = df_cur[key].max()
                        stats["mean"] = df_cur[key].mean()
                        stats["median"] = df_cur[key].median()
                        stats["mode"] = df_cur[key].mode()[0]
                        stats["variance"] = df_cur[key].var()
                        stats["count"] = df_cur[key].count()
                        stats["missing"] = df_cur[key].isna().sum()
                        stats["total"] = len(df_cur[key])
                        stats["missing_fraction"] = stats["missing"] / stats["total"]
                        df_stats = pd.concat([df_stats, pd.DataFrame([stats], index=[key])])
                
                # Save parent dataframe to CSV on disk
                fpath = os.path.join(args.output_folder, args.id,
                            f"{fpath_prefix}_{Interpretation.CONTINUOUS}_{df_str}.csv")
                df_stats.to_csv(fpath)
                logging.info(f"Saved summary stats of {Interpretation.CONTINUOUS} tmaps to {fpath}")

    # Check if any tmaps are language (strings)
    if Interpretation.LANGUAGE in [tm.interpretation for tm in tmaps]:
        for df_cur, df_str in zip([df, df.dropna()], ["union", "intersect"]):
            df_stats = pd.DataFrame()
            if df_cur.empty:
                logging.info(f"{df_str} of tensors results in empty dataframe. Skipping calculations of {Interpretation.LANGUAGE} summary statistics")
            else:
                for tm in [tm for tm in tmaps if tm.interpretation is Interpretation.LANGUAGE]:
                    if tm.channel_map:
                        for cm in tm.channel_map:
                            stats = dict()
                            key = (tm.name, cm)
                            stats["count"] = df_cur[key].count()
                            stats["count_unique"] = len(df_cur[key].value_counts())
                            stats["missing"] = df_cur[key].isna().sum()
                            stats["total"] = len(df_cur[key])
                            stats["missing_fraction"] = stats["missing"] / stats["total"]
                            df_stats = pd.concat([df_stats, pd.DataFrame([stats], index=[cm])])
                    else:
                        stats = dict()
                        key = tm.name
                        stats["count"] = df_cur[key].count()
                        stats["count_unique"] = len(df_cur[key].value_counts())
                        stats["missing"] = df_cur[key].isna().sum()
                        stats["total"] = len(df_cur[key])
                        stats["missing_fraction"] = stats["missing"] / stats["total"]
                        df_stats = pd.concat([df_stats, pd.DataFrame([stats], index=[tm.name])])

                # Save parent dataframe to CSV on disk
                fpath = os.path.join(args.output_folder, args.id,
                            f"{fpath_prefix}_{Interpretation.LANGUAGE}_{df_str}.csv")
                df_stats.to_csv(fpath)
                logging.info(f"Saved summary stats of {Interpretation.LANGUAGE} tmaps to {fpath}")


def train_multimodal_multitask(args):
    generate_train, generate_valid, generate_test = test_train_valid_tensor_generators(**args.__dict__)
    if args.variational:  # TODO: Save the encoders and decoders
        model, _, _ = make_variational_multimodal_multitask_model(**args.__dict__)
    else:
        model = make_multimodal_multitask_model(**args.__dict__)
    model = train_model_from_generators(
        model, generate_train, generate_valid, args.training_steps, args.validation_steps, args.batch_size,
        args.epochs, args.patience, args.output_folder, args.id, args.inspect_model, args.inspect_show_labels,
        anneal_rate=args.anneal_rate, anneal_shift=args.anneal_shift, anneal_max=args.anneal_max,
    )

    out_path = os.path.join(args.output_folder, args.id + '/')
    test_data, test_labels, test_paths = big_batch_from_minibatch_generator(generate_test, args.test_steps)
    return _predict_and_evaluate(model, test_data, test_labels, args.tensor_maps_in, args.tensor_maps_out, args.batch_size, args.hidden_layer, out_path, test_paths, args.alpha)


def test_multimodal_multitask(args):
    _, _, generate_test = test_train_valid_tensor_generators(**args.__dict__)
    if args.variational:
        model, _, _ = make_variational_multimodal_multitask_model(**args.__dict__)
    else:
        model = make_multimodal_multitask_model(**args.__dict__)
    out_path = os.path.join(args.output_folder, args.id + '/')
    data, labels, paths = big_batch_from_minibatch_generator(generate_test, args.test_steps)
    return _predict_and_evaluate(model, data, labels, args.tensor_maps_in, args.tensor_maps_out, args.batch_size, args.hidden_layer, out_path, paths, args.alpha)


def test_multimodal_scalar_tasks(args):
    _, _, generate_test = test_train_valid_tensor_generators(**args.__dict__)
    model = make_multimodal_multitask_model(**args.__dict__)
    p = os.path.join(args.output_folder, args.id + '/')
    return _predict_scalars_and_evaluate_from_generator(model, generate_test, args.tensor_maps_out, args.test_steps, args.hidden_layer, p, args.alpha)


def compare_multimodal_multitask_models(args):
    _, _, generate_test = test_train_valid_tensor_generators(**args.__dict__)
    models_inputs_outputs = get_model_inputs_outputs(args.model_files, args.tensor_maps_in, args.tensor_maps_out)
    input_data, output_data, paths = big_batch_from_minibatch_generator(generate_test, args.test_steps)
    common_outputs = _get_common_outputs(models_inputs_outputs, 'output')
    predictions = _get_predictions(args, models_inputs_outputs, input_data, common_outputs, 'input', 'output')
    _calculate_and_plot_prediction_stats(args, predictions, output_data, paths)


def compare_multimodal_scalar_task_models(args):
    _, _, generate_test = test_train_valid_tensor_generators(**args.__dict__)
    models_io = get_model_inputs_outputs(args.model_files, args.tensor_maps_in, args.tensor_maps_out)
    outs = _get_common_outputs(models_io, "output")
    predictions, labels, paths = _scalar_predictions_from_generator(args, models_io, generate_test, args.test_steps, outs, "input", "output")
    _calculate_and_plot_prediction_stats(args, predictions, labels, paths)


def _make_tmap_nan_on_fail(tmap):
    """
    Builds a copy TensorMap with a tensor_from_file that returns nans on errors instead of raising an error
    """
    new_tmap = copy.deepcopy(tmap)
    new_tmap.validator = lambda _, x: x  # prevent failure caused by validator

    def _tff(tm, hd5, dependents=None):
        try:
            return tmap.tensor_from_file(tm, hd5, dependents)
        except (IndexError, KeyError, ValueError, OSError, RuntimeError):
            return np.full(shape=tm.shape, fill_value=np.nan)

    new_tmap.tensor_from_file = _tff
    return new_tmap


def infer_multimodal_multitask(args):
    stats = Counter()
    tensor_paths_inferred = {}
    inference_tsv = os.path.join(args.output_folder, args.id, 'inference_' + args.id + '.tsv')
    tensor_paths = [args.tensors + tp for tp in sorted(os.listdir(args.tensors)) if os.path.splitext(tp)[-1].lower() == TENSOR_EXT]
    if args.variational:
        model, encoder, decoder = make_variational_multimodal_multitask_model(**args.__dict__)
    else:
        model = make_multimodal_multitask_model(**args.__dict__)
    no_fail_tmaps_out = [_make_tmap_nan_on_fail(tmap) for tmap in args.tensor_maps_out]
    # hard code batch size to 1 so we can iterate over file names and generated tensors together in the tensor_paths for loop
    generate_test = TensorGenerator(1, args.tensor_maps_in, no_fail_tmaps_out, tensor_paths, num_workers=0,
                                    cache_size=0, keep_paths=True, mixup=args.mixup_alpha)
    with open(inference_tsv, mode='w') as inference_file:
        # TODO: csv.DictWriter is much nicer for this
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
        inference_writer.writerow(header)

        while True:
            input_data, true_label, tensor_path = next(generate_test)
            if tensor_path[0] in tensor_paths_inferred:
                logging.info(f"Inference on {stats['count']} tensors finished. Inference TSV file at: {inference_tsv}")
                break

            prediction = model.predict(input_data)
            if len(no_fail_tmaps_out) == 1:
                prediction = [prediction]

            csv_row = [os.path.basename(tensor_path[0]).replace(TENSOR_EXT, '')]  # extract sample id
            for y, tm in zip(prediction, no_fail_tmaps_out):
                if len(tm.shape) == 1 and tm.is_continuous():
                    csv_row.append(str(tm.rescale(y)[0][0]))  # first index into batch then index into the 1x1 structure
                    if ((tm.sentinel is not None and tm.sentinel == true_label[tm.output_name()][0][0])
                            or np.isnan(true_label[tm.output_name()][0][0])):
                        csv_row.append("NA")
                    else:
                        csv_row.append(str(tm.rescale(true_label[tm.output_name()])[0][0]))
                elif len(tm.shape) == 1 and tm.is_categorical():
                    for k, i in tm.channel_map.items():
                        csv_row.append(str(y[0][tm.channel_map[k]]))
                        actual = true_label[tm.output_name()][0][i]
                        csv_row.append("NA" if np.isnan(actual) else str(actual))

            inference_writer.writerow(csv_row)
            tensor_paths_inferred[tensor_path[0]] = True
            stats['count'] += 1
            if stats['count'] % 250 == 0:
                logging.info(f"Wrote:{stats['count']} rows of inference.  Last tensor:{tensor_path[0]}")


def infer_hidden_layer_multimodal_multitask(args):
    stats = Counter()
    args.num_workers = 0
    inference_tsv = os.path.join(args.output_folder, args.id, 'hidden_inference_' + args.id + '.tsv')
    tensor_paths = [args.tensors + tp for tp in sorted(os.listdir(args.tensors)) if os.path.splitext(tp)[-1].lower() == TENSOR_EXT]
    # hard code batch size to 1 so we can iterate over file names and generated tensors together in the tensor_paths for loop
    generate_test = TensorGenerator(1, args.tensor_maps_in, args.tensor_maps_out, tensor_paths, num_workers=0,
                                    cache_size=args.cache_size, keep_paths=True, mixup=args.mixup_alpha)
    if args.variational:
        full_model, encoder, decoder = make_variational_multimodal_multitask_model(**args.__dict__)
    else:
        full_model = make_multimodal_multitask_model(**args.__dict__)
    embed_model = make_hidden_layer_model(full_model, args.tensor_maps_in, args.hidden_layer)
    dummy_input = {tm.input_name(): np.zeros((1,) + full_model.get_layer(tm.input_name()).input_shape[0][1:]) for tm in args.tensor_maps_in}
    dummy_out = embed_model.predict(dummy_input)
    latent_dimensions = int(np.prod(dummy_out.shape[1:]))
    logging.info(f'Dummy output shape is: {dummy_out.shape} latent dimensions: {latent_dimensions}')
    with open(inference_tsv, mode='w') as inference_file:
        inference_writer = csv.writer(inference_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        header = ['sample_id'] + [f'latent_{i}' for i in range(latent_dimensions)]
        inference_writer.writerow(header)

        while True:
            input_data, _, tensor_path = next(generate_test)
            if tensor_path[0] in stats:
                logging.info(f"Latent space inference on {stats['count']} tensors finished. Inference TSV file at: {inference_tsv}")
                break

            sample_id = os.path.basename(tensor_path[0]).replace(TENSOR_EXT, '')
            prediction = embed_model.predict(input_data)
            prediction = np.reshape(prediction, (latent_dimensions,))
            csv_row = [sample_id] + [f'{prediction[i]}' for i in range(latent_dimensions)]
            inference_writer.writerow(csv_row)
            stats[tensor_path[0]] += 1
            stats['count'] += 1
            if stats['count'] % 500 == 0:
                logging.info(f"Wrote:{stats['count']} rows of latent space inference.  Last tensor:{tensor_path[0]}")


def train_shallow_model(args):
    generate_train, generate_valid, generate_test = test_train_valid_tensor_generators(**args.__dict__)
    model = make_shallow_model(args.tensor_maps_in, args.tensor_maps_out, args.learning_rate, args.model_file, args.model_layers)
    model = train_model_from_generators(model, generate_train, generate_valid, args.training_steps, args.validation_steps, args.batch_size,
                                        args.epochs, args.patience, args.output_folder, args.id, args.inspect_model, args.inspect_show_labels)

    p = os.path.join(args.output_folder, args.id + '/')
    test_data, test_labels, test_paths = big_batch_from_minibatch_generator(generate_test, args.test_steps)
    return _predict_and_evaluate(model, test_data, test_labels, args.tensor_maps_in, args.tensor_maps_out, args.batch_size, args.hidden_layer, p, test_paths, args.alpha)


def train_char_model(args):
    base_model = make_multimodal_multitask_model(**args.__dict__)
    model, char_model = make_character_model_plus(args.tensor_maps_in, args.tensor_maps_out, args.learning_rate, base_model, args.model_layers)
    generate_train, generate_valid, generate_test = test_train_valid_tensor_generators(**args.__dict__)

    model = train_model_from_generators(model, generate_train, generate_valid, args.training_steps, args.validation_steps, args.batch_size,
                                        args.epochs, args.patience, args.output_folder, args.id, args.inspect_model, args.inspect_show_labels)
    test_batch, _, test_paths = next(generate_test)
    sample_from_char_model(char_model, test_batch, test_paths)

    output_path = os.path.join(args.output_folder, args.id + '/')
    data, labels, paths = big_batch_from_minibatch_generator(generate_test, args.test_steps)
    return _predict_and_evaluate(model, data, labels, args.tensor_maps_in, args.tensor_maps_out, args.batch_size, args.hidden_layer, output_path, paths, args.alpha)


def train_siamese_model(args):
    base_model = make_multimodal_multitask_model(**args.__dict__)
    siamese_model = make_siamese_model(base_model, **args.__dict__)
    generate_train, generate_valid, generate_test = test_train_valid_tensor_generators(**args.__dict__, siamese=True)
    siamese_model = train_model_from_generators(siamese_model, generate_train, generate_valid, args.training_steps, args.validation_steps, args.batch_size,
                                                args.epochs, args.patience, args.output_folder, args.id, args.inspect_model, args.inspect_show_labels)

    data, labels, paths = big_batch_from_minibatch_generator(generate_test, args.test_steps)
    prediction = siamese_model.predict(data)
    return plot_roc_per_class(prediction, labels['output_siamese'], {'random_siamese_verification_task': 0}, args.id, os.path.join(args.output_folder, args.id + '/'))


def plot_predictions(args):
    _, _, generate_test = test_train_valid_tensor_generators(**args.__dict__)
    model = make_multimodal_multitask_model(**args.__dict__)
    data, labels, paths = big_batch_from_minibatch_generator(generate_test, args.test_steps)
    predictions = model.predict(data, batch_size=args.batch_size)
    if len(args.tensor_maps_out) == 1:
        predictions = [predictions]
    folder = os.path.join(args.output_folder, args.id, 'prediction_pngs/')
    predictions_to_pngs(predictions, args.tensor_maps_in, args.tensor_maps_out, data, labels, paths, folder)


def plot_while_training(args):
    generate_train, _, generate_test = test_train_valid_tensor_generators(**args.__dict__)
    test_data, test_labels, test_paths = big_batch_from_minibatch_generator(generate_test, args.test_steps)
    model = make_multimodal_multitask_model(**args.__dict__)

    plot_folder = os.path.join(args.output_folder, args.id, 'training_frames/')
    plot_while_learning(model, args.tensor_maps_in, args.tensor_maps_out, generate_train, test_data, test_labels, test_paths, args.epochs,
                        args.batch_size, args.training_steps, plot_folder, args.write_pngs)


def saliency_maps(args):
    import tensorflow as tf
    tf.compat.v1.disable_eager_execution()
    _, _, generate_test = test_train_valid_tensor_generators(**args.__dict__)
    model = make_multimodal_multitask_model(**args.__dict__)
    data, labels, paths = big_batch_from_minibatch_generator(generate_test, args.test_steps)
    in_tensor = data[args.tensor_maps_in[0].input_name()]
    for tm in args.tensor_maps_out:
        if len(tm.shape) > 1:
            continue
        for channel in tm.channel_map:
            gradients = saliency_map(in_tensor, model, tm.output_name(), tm.channel_map[channel])
            plot_saliency_maps(in_tensor, gradients, os.path.join(args.output_folder, f'{args.id}/saliency_maps/{tm.name}_{channel}'))


def _predict_and_evaluate(model, test_data, test_labels, tensor_maps_in, tensor_maps_out, batch_size, hidden_layer, plot_path, test_paths, alpha):
    layer_names = [layer.name for layer in model.layers]
    performance_metrics = {}
    scatters = []
    rocs = []

    y_predictions = model.predict(test_data, batch_size=batch_size)
    for y, tm in zip(y_predictions, tensor_maps_out):
        if tm.output_name() not in layer_names:
            continue
        if not isinstance(y_predictions, list):  # When models have a single output model.predict returns a ndarray otherwise it returns a list
            y = y_predictions
        y_truth = np.array(test_labels[tm.output_name()])
        performance_metrics.update(evaluate_predictions(tm, y, y_truth, tm.name, plot_path, test_paths, rocs=rocs, scatters=scatters))

    if len(rocs) > 1:
        subplot_rocs(rocs, plot_path)
    if len(scatters) > 1:
        subplot_scatters(scatters, plot_path)

    test_labels_1d = {tm: np.array(test_labels[tm.output_name()]) for tm in tensor_maps_out if tm.output_name() in test_labels}
    _tsne_wrapper(model, hidden_layer, alpha, plot_path, test_paths, test_labels_1d, test_data=test_data, tensor_maps_in=tensor_maps_in, batch_size=batch_size)

    return performance_metrics


def _predict_scalars_and_evaluate_from_generator(model, test_generator, tensor_maps_out, steps, hidden_layer, plot_path, alpha):
    layer_names = [layer.name for layer in model.layers]
    model_predictions = [tm.output_name() for tm in tensor_maps_out if tm.output_name() in layer_names]
    scalar_predictions = {tm.output_name(): [] for tm in tensor_maps_out if len(tm.shape) == 1 and tm.output_name() in layer_names}
    test_labels = {tm.output_name(): [] for tm in tensor_maps_out if len(tm.shape) == 1}

    logging.info(f" in scalar predict {model_predictions} scalar predict names: {scalar_predictions.keys()} test labels: {test_labels.keys()}")
    embeddings = []
    test_paths = []
    for i in range(steps):
        batch_data, batch_labels, batch_paths = next(test_generator)
        y_predictions = model.predict(batch_data)
        test_paths.extend(batch_paths)
        if hidden_layer in layer_names:
            x_embed = embed_model_predict(model, args.tensor_maps_in, hidden_layer, batch_data, 2)
            embeddings.extend(np.copy(np.reshape(x_embed, (x_embed.shape[0], np.prod(x_embed.shape[1:])))))

        for tm_output_name in test_labels:
            test_labels[tm_output_name].extend(np.copy(batch_labels[tm_output_name]))

        for y, tm_output_name in zip(y_predictions, model_predictions):
            if not isinstance(y_predictions, list):  # When models have a single output model.predict returns a ndarray otherwise it returns a list
                y = y_predictions
            if tm_output_name in scalar_predictions:
                scalar_predictions[tm_output_name].extend(np.copy(y))

    performance_metrics = {}
    scatters = []
    rocs = []
    for tm in tensor_maps_out:
        if tm.output_name() in scalar_predictions:
            y_predict = np.array(scalar_predictions[tm.output_name()])
            y_truth = np.array(test_labels[tm.output_name()])
            performance_metrics.update(evaluate_predictions(tm, y_predict, y_truth, tm.name, plot_path, test_paths, rocs=rocs, scatters=scatters))

    if len(rocs) > 1:
        subplot_rocs(rocs, plot_path)
    if len(scatters) > 1:
        subplot_scatters(scatters, plot_path)
    if len(embeddings) > 0:
        test_labels_1d = {tm: np.array(test_labels[tm.output_name()]) for tm in tensor_maps_out if tm.output_name() in test_labels}
        _tsne_wrapper(model, hidden_layer, alpha, plot_path, test_paths, test_labels_1d, embeddings=embeddings)

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
        args.tensor_maps_out = models_inputs_outputs[model_file][output_prefix]
        args.tensor_maps_in = models_inputs_outputs[model_file][input_prefix]
        args.model_file = model_file
        model = make_multimodal_multitask_model(**args.__dict__)
        model_name = os.path.basename(model_file).replace(MODEL_EXT, '_')

        # We can feed 'model.predict()' the entire input data because it knows what subset to use
        y_pred = model.predict(input_data, batch_size=args.batch_size)

        for i, tm in enumerate(args.tensor_maps_out):
            if tm in outputs:
                if len(args.tensor_maps_out) == 1:
                    predictions[tm][model_name] = y_pred
                else:
                    predictions[tm][model_name] = y_pred[i]

    return predictions


def _scalar_predictions_from_generator(args, models_inputs_outputs, generator, steps, outputs, input_prefix, output_prefix):
    """Makes multi-modal scalar predictions for a given number of models.

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
    models = {}
    test_paths = []
    scalar_predictions = {}
    test_labels = {tm.output_name(): [] for tm in args.tensor_maps_out if len(tm.shape) == 1}

    for model_file in models_inputs_outputs:
        args.model_file = model_file
        args.tensor_maps_in = models_inputs_outputs[model_file][input_prefix]
        args.tensor_maps_out = models_inputs_outputs[model_file][output_prefix]
        model = make_multimodal_multitask_model(**args.__dict__)
        model_name = os.path.basename(model_file).replace(MODEL_EXT, '')
        models[model_name] = model
        scalar_predictions[model_name] = [tm for tm in models_inputs_outputs[model_file][output_prefix] if len(tm.shape) == 1]

    predictions = defaultdict(dict)
    for j in range(steps):
        input_data, labels, paths = next(generator)
        test_paths.extend(paths)
        for tl in test_labels:
            test_labels[tl].extend(np.copy(labels[tl]))

        for model_name, model_file in zip(models, models_inputs_outputs):
            # We can feed 'model.predict()' the entire input data because it knows what subset to use
            y_predictions = models[model_name].predict(input_data)

            for y, tm in zip(y_predictions, models_inputs_outputs[model_file][output_prefix]):
                if not isinstance(y_predictions, list):  # When models have a single output model.predict returns a ndarray otherwise it returns a list
                    y = y_predictions
                if j == 0 and tm in scalar_predictions[model_name]:
                    predictions[tm][model_name] = []
                if tm in scalar_predictions[model_name]:
                    predictions[tm][model_name].extend(np.copy(y))

    for tm in predictions:
        logging.info(f"{tm.output_name()} labels: {len(test_labels[tm.output_name()])}")
        test_labels[tm.output_name()] = np.array(test_labels[tm.output_name()])
        for m in predictions[tm]:
            logging.info(f"{tm.output_name()} model: {m} prediction length:{len(predictions[tm][m])}")
            predictions[tm][m] = np.array(predictions[tm][m])

    return predictions, test_labels, test_paths


def _calculate_and_plot_prediction_stats(args, predictions, outputs, paths):
    rocs = []
    scatters = []
    for tm in args.tensor_maps_out:
        if tm not in predictions:
            continue
        plot_title = tm.name+'_'+args.id
        plot_folder = os.path.join(args.output_folder, args.id)

        if tm.is_categorical() and tm.axes() == 1:
            msg = "For tm '{}' with channel map {}: sum truth = {}; sum pred = {}"
            for m in predictions[tm]:
                logging.info(msg.format(tm.name, tm.channel_map, np.sum(outputs[tm.output_name()], axis=0), np.sum(predictions[tm][m], axis=0)))
            plot_rocs(predictions[tm], outputs[tm.output_name()], tm.channel_map, plot_title, plot_folder)
            rocs.append((predictions[tm], outputs[tm.output_name()], tm.channel_map))
        elif tm.is_categorical() and tm.axes() == 4:
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
        elif tm.is_continuous() and tm.axes() == 1:
            scaled_predictions = {k: tm.rescale(predictions[tm][k]) for k in predictions[tm]}
            plot_scatters(scaled_predictions, tm.rescale(outputs[tm.output_name()]), plot_title, plot_folder, paths)
            scatters.append((scaled_predictions, tm.rescale(outputs[tm.output_name()]), plot_title, None))
            coefs = get_pearson_coefficients(scaled_predictions, tm.rescale(outputs[tm.output_name()]))
            log_pearson_coefficients(coefs, tm.name)
        else:
            scaled_predictions = {k: tm.rescale(predictions[tm][k]) for k in predictions[tm]}
            plot_scatters(scaled_predictions, tm.rescale(outputs[tm.output_name()]), plot_title, plot_folder)
            coefs = get_pearson_coefficients(scaled_predictions, tm.rescale(outputs[tm.output_name()]))
            log_pearson_coefficients(coefs, tm.name)

    if len(rocs) > 1:
        subplot_comparison_rocs(rocs, plot_folder)
    if len(scatters) > 1:
        subplot_comparison_scatters(scatters, plot_folder)


def _tsne_wrapper(model, hidden_layer_name, alpha, plot_path, test_paths, test_labels, test_data=None, tensor_maps_in=None, batch_size=16, embeddings=None):
    """Plot 2D t-SNE of a model's hidden layer colored by many different co-variates.

    Callers must provide either model's embeddings or test_data on which embeddings will be inferred

    :param model: Keras model
    :param hidden_layer_name: String name of the hidden layer whose embeddings will be visualized
    :param alpha: Transparency of each data point
    :param plot_path: Image file name and path for the t_SNE plot
    :param test_paths: Paths for hd5 file containing each sample
    :param test_labels: Dictionary mapping TensorMaps to numpy arrays of labels (co-variates) to color code the t-SNE plots with
    :param test_data: Input data for the model necessary if embeddings is None
    :param tensor_maps_in: Input TensorMaps of the model necessary if embeddings is None
    :param batch_size: Batch size necessary if embeddings is None
    :param embeddings: (optional) Model's embeddings
    :return: None
    """
    if hidden_layer_name not in [layer.name for layer in model.layers]:
        logging.warning(f"Can't compute t-SNE, layer:{hidden_layer_name} not in provided model.")
        return

    if embeddings is None:
        embeddings = embed_model_predict(model, tensor_maps_in, hidden_layer_name, test_data, batch_size)

    gene_labels = []
    label_dict, categorical_labels, continuous_labels = test_labels_to_label_map(test_labels, len(test_paths))
    if len(categorical_labels) > 0 or len(continuous_labels) > 0 or len(gene_labels) > 0:
        plot_tsne(embeddings, categorical_labels, continuous_labels, gene_labels, label_dict, plot_path, alpha)


if __name__ == '__main__':
    args = parse_args()
    run(args)  # back to the top
