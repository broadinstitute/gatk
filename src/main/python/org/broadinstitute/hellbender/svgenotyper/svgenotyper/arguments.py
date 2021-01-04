import argparse

"""
Parse command line arguments.
"""


def parse_args_train():
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', help='Input VCF path', required=True)
    parser.add_argument('--coverage-file', help='Table of sample mean per-base coverage', required=True)
    parser.add_argument('--samples-file', help='List of sample ids', required=True)
    parser.add_argument('--output-name', help='Output VCF base name', required=True)
    parser.add_argument('--output-dir', help='Output VCF directory', default='.')

    parser.add_argument('--device', default='cpu', help='Device for Torch backend (e.g. "cpu", "cuda")')
    parser.add_argument('--random-seed', default=92837488, type=int, help='PRNG seed')

    parser.add_argument('--num-states', default=None, type=int, help='Max number of genotype states')
    parser.add_argument('--depth-dilution-factor', default=1e-6, type=float, help='Dilution factor for depth posteriors')
    parser.add_argument('--eps-pe', default=0.1, type=float, help='Mean of PE noise prior')
    parser.add_argument('--eps-sr1', default=0.1, type=float, help='Mean of SR1 noise prior')
    parser.add_argument('--eps-sr2', default=0.1, type=float, help='Mean of SR2 noise prior')
    parser.add_argument('--lambda-pe', default=0.1, type=float, help='Mean of PE variance factor prior')
    parser.add_argument('--lambda-sr1', default=0.1, type=float, help='Mean of SR1 variance factor prior')
    parser.add_argument('--lambda-sr2', default=0.1, type=float, help='Mean of SR2 variance factor prior')
    parser.add_argument('--phi-pe', default=0.1, type=float, help='Variance of PE bias prior')
    parser.add_argument('--phi-sr1', default=0.1, type=float, help='Variance of SR1 bias prior')
    parser.add_argument('--phi-sr2', default=0.1, type=float, help='Variance of SR2 bias prior')
    parser.add_argument('--eta-q', default=0.1, type=float, help='Mean of allele frequency prior')
    parser.add_argument('--eta-r', default=0.01, type=float, help='Mean of secondary allele prior (DUP only)')

    parser.add_argument('--lr-decay', default=1000., type=float, help='Learning rate decay constant (lower is faster)')
    parser.add_argument('--lr-min', default=1e-3, type=float, help='Minimum learning rate')
    parser.add_argument('--lr-init', default=0.01, type=float, help='Initial learning rate')
    parser.add_argument('--adam-beta1', default=0.9, type=float, help='ADAM beta1 constant')
    parser.add_argument('--adam-beta2', default=0.999, type=float, help='ADAM beta2 constant')

    parser.add_argument('--max-iter', default=2000, type=int, help='Max number of training iterations')
    parser.add_argument('--iter-log-freq', default=100, type=int, help='Number of iterations between log messages')
    parser.add_argument('--jit', action='store_true', help='Enable JIT compilation')

    args = parser.parse_args()

    print('Arguments are', args)
    return args


def parse_args_genotype():
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', help='Input VCF path', required=True)
    parser.add_argument('--model-dir', help='Input model directory', required=True)
    parser.add_argument('--model-name', help='Input model name', required=True)
    parser.add_argument('--output', help='Output VCF path', required=True)
    parser.add_argument('--plots-dir', default=None, help='If specified, generates plots in this directory')

    parser.add_argument('--device', default='cpu', help='Device for Torch backend (e.g. "cpu", "cuda")')
    parser.add_argument('--random-seed', default=92837488, type=int, help='PRNG seed')

    parser.add_argument('--jit', action='store_true', help='Enable JIT compilation')

    parser.add_argument('--genotype-predictive-samples', default=1000, type=int,
                        help='Number of samples for predictive inference')
    parser.add_argument('--genotype-discrete-samples', default=1000, type=int,
                        help='Number of samples for discrete inference')
    parser.add_argument('--genotype-discrete-log-freq', default=100, type=int,
                        help='Number of discrete inference samples between log messages')

    args = parser.parse_args()

    print('Arguments are', args)
    return args
