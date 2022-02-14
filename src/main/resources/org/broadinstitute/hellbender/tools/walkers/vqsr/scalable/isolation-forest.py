import argparse
import h5py
import sklearn.ensemble
import sklearn.impute
import numpy as np
import dill


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


def do_work(annotations_file,
            hyperparameters_json_file,
            output_prefix):
    print('Reading...')
    annotation_names_i, X_ni = read_annotations(annotations_file)

    print('Imputing...')
    imputer = sklearn.impute.SimpleImputer(strategy='median')
    imputed_X_ni = imputer.fit_transform(X_ni)

    # TODO ingest hyperparameters

    print(f'Training IsolationForest with {imputed_X_ni.shape[0]} training sites x {imputed_X_ni.shape[1]} annotations...')
    clf = sklearn.ensemble.IsolationForest(random_state=0, verbose=True)
    clf.fit(imputed_X_ni)

    scorer_lambda = lambda test_X_ni: clf.score_samples(imputer.transform(test_X_ni))

    print(f'Pickling scorer...')
    output_scorer_pkl_file = f'{output_prefix}.scorer.pkl'
    with open(output_scorer_pkl_file, 'wb') as f:
        dill.dump(scorer_lambda, f)
    print(f'Scorer pickled to {output_scorer_pkl_file}.')

    print(f'Pickling imputer...')
    output_imputer_pkl_file = f'{output_prefix}.imputer.pkl'
    with open(output_imputer_pkl_file, 'wb') as f:
        dill.dump(imputer, f)
    print(f'Imputer pickled to {output_imputer_pkl_file}.')

    print(f'Pickling model...')
    output_model_pkl_file = f'{output_prefix}.model.pkl'
    with open(output_model_pkl_file, 'wb') as f:
        dill.dump(clf, f)
    print(f'Model pickled to {output_model_pkl_file}.')

    # print('Scoring...')
    # scores_n = clf.score_samples(imputed_X_ni)
    #
    # print(f'Writing scores...')
    # output_training_scores_file = f'{output_prefix}.trainingScores.hdf5'
    # with h5py.File(output_training_scores_file, 'w') as f:
    #     scores_dset = f.create_dataset('data/scores', (len(scores_n[is_training_n]),), dtype='d')
    #     scores_dset[:] = scores_n[is_training_n]
    # print(f'Training scores written to {output_training_scores_file}.')
    # output_truth_scores_file = f'{output_prefix}.truthScores.hdf5'
    # with h5py.File(output_truth_scores_file, 'w') as f:
    #     scores_dset = f.create_dataset('data/scores', (len(scores_n[is_truth_n]),), dtype='d')
    #     scores_dset[:] = scores_n[is_truth_n]
    # print(f'Truth scores written to {output_truth_scores_file}.')


def main():
    parser = argparse.ArgumentParser()

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

    annotations_file = args.annotations_file
    hyperparameters_json_file = args.hyperparameters_json_file
    output_prefix = args.output_prefix

    do_work(annotations_file,
            hyperparameters_json_file,
            output_prefix)


if __name__ == '__main__':
    main()
