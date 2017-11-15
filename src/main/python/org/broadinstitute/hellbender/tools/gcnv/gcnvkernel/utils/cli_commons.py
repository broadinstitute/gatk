import argparse
import logging
import sys

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


def set_logging_config_from_args(parsed_args):
    """Configures python logger according to parsed arguments."""
    logging.basicConfig(level=log_level_map[parsed_args.logfile_log_level],
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=parsed_args.logfile if hasattr(parsed_args, 'logfile') else '/dev/null',
                        filemode='w')
    console = logging.StreamHandler(stream=sys.stdout)
    console.setLevel(log_level_map[parsed_args.console_log_level])
    formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
