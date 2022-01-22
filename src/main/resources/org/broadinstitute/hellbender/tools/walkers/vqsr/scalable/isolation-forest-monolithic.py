import argparse
import h5py
import sklearn.ensemble
import sklearn.impute
import numpy as np


def read_annotations(h5file):
    with h5py.File(h5file, 'r') as f:
        annotation_names_i = f['/data/annotation_names'][:].astype(str)
        X_ni = f['/data/annotations'][:]
        is_anti_training_n = f['/data/is_anti_training'][:].astype(bool)
        is_training_n = f['/data/is_training'][:].astype(bool)
        is_truth_n = f['/data/is_truth'][:].astype(bool)
    return annotation_names_i, X_ni, is_anti_training_n, is_training_n, is_truth_n


def do_work(raw_annotations_file,
            annotations_file,
            hyperparameters_json_file,
            output_prefix):
    annotation_names_i, X_ni, is_anti_training_n, is_training_n, is_truth_n = read_annotations(raw_annotations_file)

    imputer = sklearn.impute.SimpleImputer(strategy='median')
    imputer.fit(X_ni[is_training_n])
    imputed_X_ni = imputer.transform(X_ni)
    training_X_ni = imputed_X_ni[is_training_n]

    print('Number of training sites =', np.sum(is_training_n))

    # TODO plot

    # TODO ingest hyperparameters

    print('Training IsolationForest...')
    clf = sklearn.ensemble.IsolationForest(random_state=0, verbose=True)
    clf.fit(training_X_ni)

    print('Scoring...')
    scores = clf.score_samples(imputed_X_ni)

    outputScoresFile = f'{output_prefix}.scores.hdf5'
    with h5py.File(outputScoresFile, 'w') as f:
        scores_dset = f.create_dataset('scores', (len(scores),), dtype='d')
        scores_dset[:] = scores


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--raw_annotations_file',
                        type=str,
                        required=True,
                        help='')

    parser.add_argument('--annotations_file',
                        type=str,
                        required=True,
                        help='')

    parser.add_argument('--hyperparameters_json_file',
                        type=str,
                        required=True,
                        help='')

    parser.add_argument('--output_prefix',
                        type=str,
                        required=True,
                        help='')

    args = parser.parse_args()

    raw_annotations_file = args.raw_annotations_file
    annotations_file = args.annotations_file
    hyperparameters_json_file = args.hyperparameters_json_file
    output_prefix = args.output_prefix

    do_work(raw_annotations_file,
            annotations_file,
            hyperparameters_json_file,
            output_prefix)


if __name__ == '__main__':
    main()
