import os
import argparse
import logging
import sys
import importlib
import datetime as dt

logging_level = logging.DEBUG
console_logging_level = logging.INFO


class GCNVHelpFormatter(argparse.HelpFormatter):
    """Argument parser help formatter for gcnvkernel CLI scripts."""

    def _get_help_string(self, action):
        _help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    _help += ' (default: %(default)s)'
        return _help

    def _get_default_metavar_for_optional(self, action):
        return action.type.__name__

    def _get_default_metavar_for_positional(self, action):
        return action.type.__name__


class GATKStyleLogFormatter(logging.Formatter):
    def formatTime(self, record, datefmt=None):
        ct = dt.datetime.fromtimestamp(record.created)
        t = ct.strftime("%H:%M:%S")
        s = "%s.%03d" % (t, record.msecs)
        return s


def set_logging_config():
    """Configures python logger."""
    logging.basicConfig(level=logging_level,
                        format='%(asctime)s %(levelname)s %(name)s - %(message)s',
                        filename='/dev/null',
                        filemode='w')
    console = logging.StreamHandler(stream=sys.stdout)
    console.setLevel(console_logging_level)
    console.setFormatter(GATKStyleLogFormatter(
        '%(asctime)s %(levelname)s %(name)s - %(message)s', datefmt=None))
    logging.getLogger('').addHandler(console)


def set_theano_flags(parser, logger):
    """Set THEANO_FLAGS environment variable. Rightmost flags take precedence, so users can override
    or specify additional flags by setting THEANO_FLAGS before running gCNV scripts.  The argument
    basecompile_dir is added to the parser as a side effect.  For details, see documentation at
    http://deeplearning.net/software/theano/library/config.html."""

    parser.add_argument("--base_compiledir",
                        type=str,
                        required=False,
                        default=None,
                        help="Specifies the base directory in which theano compilation subdirectories will be placed.")

    args, unknown_args = parser.parse_known_args()

    # set flags
    user_theano_flags = os.environ.get("THEANO_FLAGS")
    default_theano_flags = "device=cpu,floatX=float64,optimizer=fast_run,compute_test_value=ignore," + \
                           "openmp=true,blas.ldflags=-lmkl_rt,openmp_elemwise_minsize=10" + \
                           ("" if args.base_compiledir is None  # will set base_compiledir = $HOME/.theano
                            else ",base_compiledir={base_compiledir}".format(base_compiledir=args.base_compiledir))
    theano_flags = default_theano_flags + \
                   ("" if user_theano_flags is None else "," + user_theano_flags)
    os.environ["THEANO_FLAGS"] = theano_flags
    logger.info("THEANO_FLAGS environment variable has been set to: {theano_flags}".format(theano_flags=theano_flags))

    # reload theano to reset config
    importlib.import_module("theano")
    logger.debug(theano.config)

def str_to_bool(value: str):
    if value.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif value.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')