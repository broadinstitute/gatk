import argparse
import h5py
import sklearn.ensemble
import sklearn.impute
import numpy as np
import dill
import json
import logging


logger = logging.getLogger("isolation_forest")
logging.basicConfig(level=logging.INFO)


def read_annotations(h5file):
    with h5py.File(h5file, 'r') as f:
        annotation_names_i = f['/annotations/names'][()].astype(str)

        # read chunked annotations
        num_chunks = int(f['/annotations/num_chunks'][()])
        num_columns = int(f['/annotations/num_columns'][()])
        num_rows = int(f['/annotations/num_rows'][()])
        X_ni = np.zeros((num_rows, num_columns))
        n = 0
        for chunk_index in range(num_chunks):
            chunk_ni = f[f'/annotations/chunk_{chunk_index}'][()]
            num_rows_in_chunk = len(chunk_ni)
            X_ni[n:n + num_rows_in_chunk, :] = chunk_ni
            n += num_rows_in_chunk
        assert n == num_rows
    return annotation_names_i, X_ni


def train(annotations_file,
          hyperparameters_json_file,
          output_prefix):
    logger.info('Reading annotations...')
    annotation_names_i, X_ni = read_annotations(annotations_file)
    logger.info(f'Annotations: {annotation_names_i}.')

    logger.info('Reading hyperparameters...')
    with open(hyperparameters_json_file) as json_file:
        hyperparameters_kwargs = json.load(json_file)
    logger.info('Hyperparameters:', hyperparameters_kwargs)

    logger.info('Imputing annotations...')
    imputer = sklearn.impute.SimpleImputer(strategy='median')
    imputed_X_ni = imputer.fit_transform(X_ni)

    # SimpleImputer will drop any features that are completely missing, resulting in different shapes for
    # imputed_X_ni and X_ni and misalignment of features when training and scoring downstream if not checked.
    # We externally check for and fail in the presence of any entirely missing features, but we do a redundant check here.
    assert imputed_X_ni.shape == X_ni.shape, \
        f'Shape of imputed annotations differs from shape of raw annotations; at least one feature is completely missing ' \
        f'and hence dropped during imputation.'

    logger.info(f'Training IsolationForest with {imputed_X_ni.shape[0]} training sites x {imputed_X_ni.shape[1]} annotations...')
    clf = sklearn.ensemble.IsolationForest(**hyperparameters_kwargs)
    clf.fit(imputed_X_ni)
    logger.info('Training complete.')

    def score_samples(test_annotation_names_i,
                      test_X_ni):
        assert np.array_equal(test_annotation_names_i, annotation_names_i), \
            f'Input annotation names ({test_annotation_names_i}) must be identical to those used to train the scorer ({annotation_names_i}).'
        return clf.score_samples(imputer.transform(test_X_ni)) # TODO sklearn's implementation is single-threaded, but this could perhaps be parallelized

    scorer_lambda = lambda test_annotation_names_i, test_X_ni: score_samples(test_annotation_names_i, test_X_ni)

    logger.info(f'Pickling scorer...')
    output_scorer_pkl_file = f'{output_prefix}.scorer.pkl'
    with open(output_scorer_pkl_file, 'wb') as f:
        dill.dump(scorer_lambda, f) # the dill package can be used to pickle lambda functions
    logger.info(f'Scorer pickled to {output_scorer_pkl_file}.')


def score(annotations_file,
          scorer_pkl_file,
          output_scores_file):
    annotation_names_i, X_ni = read_annotations(annotations_file)

    with open(scorer_pkl_file, 'rb') as f:
        scorer_lambda = dill.load(f)
    score_n = scorer_lambda(annotation_names_i, X_ni)

    with h5py.File(output_scores_file, 'w') as f:
        scores_dset = f.create_dataset('data/scores', (len(score_n),), dtype='d')
        scores_dset[:] = score_n


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--annotations_file',
                        type=str,
                        required=True,
                        help='')

    parser.add_argument('--hyperparameters_json_file',
                        type=str,
                        required=False,
                        help='')

    parser.add_argument('--output_prefix',
                        type=str,
                        required=False,
                        help='')

    parser.add_argument('--scorer_pkl_file',
                        type=str,
                        required=False,
                        help='')

    parser.add_argument('--output_scores_file',
                        type=str,
                        required=False,
                        help='')

    args = parser.parse_args()

    annotations_file = args.annotations_file

    # this script can handle both training and scoring; we check the passed arguments to determine which is appropriate
    if args.hyperparameters_json_file is not None and args.output_prefix is not None and \
        args.scorer_pkl_file is None and args.output_scores_file is None:
        hyperparameters_json_file = args.hyperparameters_json_file
        output_prefix = args.output_prefix
        train(annotations_file,
              hyperparameters_json_file,
              output_prefix)
    elif args.hyperparameters_json_file is None and args.output_prefix is None and \
            args.scorer_pkl_file is not None and args.output_scores_file is not None:
        scorer_pkl_file = args.scorer_pkl_file
        output_scores_file = args.output_scores_file
        score(annotations_file,
              scorer_pkl_file,
              output_scores_file)
    else:
        raise


if __name__ == '__main__':
    main()
