import argparse
from svgenotyper import train

parser = argparse.ArgumentParser()

parser.add_argument("--variants_file", type=str, required=True)
parser.add_argument("--coverage_file", type=str, required=True)
parser.add_argument("--samples_file", type=str, required=True)
parser.add_argument("--output_name", type=str, required=True)
parser.add_argument("--output_dir", type=str, required=True)
parser.add_argument("--svtype", type=str, required=True)
parser.add_argument("--device", type=str, required=True)

parser.add_argument("--num_states", type=int, required=True)
parser.add_argument("--random_seed", type=int, required=True)
parser.add_argument("--read_length", type=int, required=True)
parser.add_argument("--max_iter", type=int, required=True)
parser.add_argument("--iter_log_freq", type=int, required=True)

parser.add_argument("--depth_dilution_factor", type=float, required=True)
parser.add_argument("--eps_pe", type=float, required=True)
parser.add_argument("--eps_sr1", type=float, required=True)
parser.add_argument("--eps_sr2", type=float, required=True)
parser.add_argument("--phi_pe", type=float, required=True)
parser.add_argument("--phi_sr1", type=float, required=True)
parser.add_argument("--phi_sr2", type=float, required=True)
parser.add_argument("--lr_decay", type=float, required=True)
parser.add_argument("--lr_min", type=float, required=True)
parser.add_argument("--lr_init", type=float, required=True)
parser.add_argument("--adam_beta1", type=float, required=True)
parser.add_argument("--adam_beta2", type=float, required=True)

parser.add_argument("--jit", action="store_true")

if __name__ == "__main__":
    # parse arguments
    args = parser.parse_args()
    #run genotyper
    train.run(args=args.__dict__)
