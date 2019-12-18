import argparse
from svgenotyper import cnv_train

parser = argparse.ArgumentParser()

parser.add_argument("--data_file", type=str, required=True)
parser.add_argument("--sample_depth_file", type=str, required=True)
parser.add_argument("--samples_file", type=str, required=True)
parser.add_argument("--output_name", type=str, required=True)
parser.add_argument("--output_dir", type=str, required=True)
parser.add_argument("--device", type=str, required=True)

parser.add_argument("--num_states", type=int, required=True)
parser.add_argument("--random_seed", type=int, required=True)
parser.add_argument("--read_length", type=int, required=True)
parser.add_argument("--max_iter", type=int, required=True)
parser.add_argument("--iter_log_freq", type=int, required=True)

parser.add_argument("--var_phi_bin", type=float, required=True)
parser.add_argument("--var_phi_sample", type=float, required=True)
parser.add_argument("--alpha_ref", type=float, required=True)
parser.add_argument("--alpha_non_ref", type=float, required=True)
parser.add_argument("--mu_eps", type=float, required=True)

parser.add_argument("--lr_decay", type=float, required=True)
parser.add_argument("--lr_min", type=float, required=True)
parser.add_argument("--lr_init", type=float, required=True)
parser.add_argument("--adam_beta1", type=float, required=True)
parser.add_argument("--adam_beta2", type=float, required=True)

parser.add_argument("--jit", action="store_true")

if __name__ == "__main__":
    args = parser.parse_args()
    cnv_train.run(args=args.__dict__)
