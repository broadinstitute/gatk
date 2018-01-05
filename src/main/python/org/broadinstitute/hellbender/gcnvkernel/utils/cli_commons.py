import argparse
import logging
import sys
import datetime as dt

log_level_map = {
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "DEBUG": logging.DEBUG
}


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


def add_logging_args_to_argparse(parser: argparse.ArgumentParser) -> None:
    """Adds logging-related arguments to a given `ArgumentParser`."""
    parser.add_argument("--console_log_level",
                        type=str,
                        choices=["INFO", "WARNING", "DEBUG"],
                        default="INFO",
                        help="Console logging verbosity level")

    parser.add_argument("--logfile_log_level",
                        type=str,
                        choices=["INFO", "WARNING", "DEBUG"],
                        default="DEBUG",
                        help="Logfile logging verbosity level")

    parser.add_argument("--logfile",
                        type=str,
                        required=False,
                        default=argparse.SUPPRESS,
                        help="If provided, the output log will be written to file as well")


class GATKStyleLogFormatter(logging.Formatter):
    def formatTime(self, record, datefmt=None):
        ct = dt.datetime.fromtimestamp(record.created)
        t = ct.strftime("%H:%M:%S")
        s = "%s.%03d" % (t, record.msecs)
        return s


def set_logging_config_from_args(parsed_args):
    """Configures python logger according to parsed arguments."""
    logging.basicConfig(level=log_level_map[parsed_args.logfile_log_level],
                        format='%(asctime)s %(levelname)s %(name)s - %(message)s',
                        filename=parsed_args.logfile if hasattr(parsed_args, 'logfile') else '/dev/null',
                        filemode='w')
    console = logging.StreamHandler(stream=sys.stdout)
    console.setLevel(log_level_map[parsed_args.console_log_level])
    console.setFormatter(GATKStyleLogFormatter(
        '%(asctime)s %(levelname)s %(name)s - %(message)s', datefmt=None))
    logging.getLogger('').addHandler(console)
