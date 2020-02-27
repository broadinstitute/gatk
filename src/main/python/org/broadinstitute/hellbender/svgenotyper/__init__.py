from .svgenotyper import train, infer
from .svgenotyper.arguments import parse_args_infer, parse_args_train
from .svgenotyper.io import write_vcf
from .svgenotyper.plot import plot_sites, plot_sites_per_genome_summary
from ._version import __version__