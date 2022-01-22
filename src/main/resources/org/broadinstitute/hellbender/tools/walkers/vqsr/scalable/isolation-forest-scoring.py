import argparse
import h5py
import numpy as np
import dill


def read_annotations(h5file):
    with h5py.File(h5file, 'r') as f:
        annotation_names_i = f['/data/annotation_names'][()].astype(str)
        is_training_n = f['/data/is_training'][()].astype(bool)
        is_truth_n = f['/data/is_truth'][()].astype(bool)

        # read chunked annotations
        num_chunks = int(f['/data/annotations/num_chunks'][()])
        num_columns = int(f['/data/annotations/num_columns'][()])
        num_rows = int(f['/data/annotations/num_rows'][()])
        X_ni = np.zeros((num_rows, num_columns))
        n = 0
        for chunk_index in range(num_chunks):
            chunk_ni = f[f'/data/annotations/chunk_{chunk_index}'][()]
            num_rows_in_chunk = len(chunk_ni)
            X_ni[n:n + num_rows_in_chunk, :] = chunk_ni
            n += num_rows_in_chunk
        assert n == num_rows
    return annotation_names_i, X_ni, is_training_n, is_truth_n


def do_work(raw_annotations_file,
            scorer_pkl_file,
            output_scores_file):
    annotation_names_i, X_ni, is_training_n, is_truth_n = read_annotations(raw_annotations_file)

    with open(scorer_pkl_file, 'rb') as f:
        scorer_lambda = dill.load(f)
    scores = scorer_lambda(X_ni)

    with h5py.File(output_scores_file, 'w') as f:
        scores_dset = f.create_dataset('scores', (len(scores),), dtype='d')
        scores_dset[:] = scores


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--raw_annotations_file',
                        type=str,
                        required=True,
                        help='')

    parser.add_argument('--scorer_pkl_file',
                        type=str,
                        required=True,
                        help='')

    parser.add_argument('--output_scores_file',
                        type=str,
                        required=True,
                        help='')

    args = parser.parse_args()

    raw_annotations_file = args.raw_annotations_file
    scorer_pkl_file = args.scorer_pkl_file
    output_scores_file = args.output_scores_file

    do_work(raw_annotations_file,
            scorer_pkl_file,
            output_scores_file)


if __name__ == '__main__':
    main()
