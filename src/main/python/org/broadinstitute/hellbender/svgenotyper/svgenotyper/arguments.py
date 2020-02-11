import argparse


def parse_args():
    """Parse command line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', help='Input VCF path')
    parser.add_argument('--coverage-file', help='Table of sample mean per-base coverage')
    parser.add_argument('--output', help='Output VCF path')
    parser.add_argument('--plots-dir', default=None, help='If specified, generates plots in this directory')

    parser.add_argument('--device', default='cpu', help='Device for Torch backend (e.g. "cpu", "cuda")')
    parser.add_argument('--random-seed', default=92837488, type=int, help='PRNG seed')

    parser.add_argument('--num-states', default=5, type=int, help='Max number of genotype copy states')
    parser.add_argument('--depth-dilution-factor', default=1e-6, type=float, help='Dilution factor for depth posteriors')
    parser.add_argument('--eps-pe', default=0.1, type=float, help='Mean of PE noise prior')
    parser.add_argument('--eps-sr1', default=0.1, type=float, help='Mean of SR1 noise prior')
    parser.add_argument('--eps-sr2', default=0.1, type=float, help='Mean of SR2 noise prior')

    parser.add_argument('--lr-decay', default=1000., type=float, help='Learning rate decay constant (lower is faster)')
    parser.add_argument('--lr-min', default=1e-4, type=float, help='Minimum learning rate')
    parser.add_argument('--lr-init', default=0.1, type=float, help='Initial learning rate')
    parser.add_argument('--adam-beta1', default=0.9, type=float, help='ADAM beta1 constant')
    parser.add_argument('--adam-beta2', default=0.999, type=float, help='ADAM beta2 constant')

    parser.add_argument('--max-iter', default=5000, type=int, help='Max number of training iterations')
    parser.add_argument('--iter-log-freq', default=100, type=int, help='Number of iterations between log messages')
    parser.add_argument('--jit', action='store_true', help='Enable JIT compilation')

    parser.add_argument('--num-infer-samples', default=100, type=int, help='Number of samples for inference')

    # Parse, print, set annotations and seed
    args = parser.parse_args()

    print('Arguments are', args)
    return args

