import argparse
from svgenotyper import cnv_infer

parser = argparse.ArgumentParser()

parser.add_argument("--output", type=str, required=True)
parser.add_argument("--model_name", type=str, required=True)
parser.add_argument("--model_dir", type=str, required=True)
parser.add_argument("--device", type=str, required=True)

parser.add_argument("--random_seed", type=int, required=True)
parser.add_argument("--infer_predictive_samples", type=int, required=True)
parser.add_argument("--infer_predictive_iter", type=int, required=True)
parser.add_argument("--infer_discrete_samples", type=int, required=True)
parser.add_argument("--infer_discrete_log_freq", type=int, required=True)

parser.add_argument("--jit", action="store_true")

if __name__ == "__main__":
    args = parser.parse_args()
    cnv_infer.run(args=args.__dict__)
