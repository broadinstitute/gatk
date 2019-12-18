import argparse
from svgenotyper import genotype

parser = argparse.ArgumentParser()

parser.add_argument("--output", type=str, required=True)
parser.add_argument("--model_name", type=str, required=True)
parser.add_argument("--model_dir", type=str, required=True)
parser.add_argument("--svtype", type=str, required=True)
parser.add_argument("--device", type=str, required=True)

parser.add_argument("--random_seed", type=int, required=True)
parser.add_argument("--genotype_discrete_samples", type=int, required=True)
parser.add_argument("--genotype_discrete_log_freq", type=int, required=True)

parser.add_argument("--jit", action="store_true")

if __name__ == "__main__":
    # parse arguments
    args = parser.parse_args()
    #run genotyper
    genotype.run(args=args.__dict__)
