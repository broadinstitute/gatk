"""Provides config settings for the logger and a way to load them"""

import sys
import os
import errno
import logging


def load_config(log_level, log_dir, log_file_basename, log_file_suffix):
    from logging import config as logging_config

    try:
        os.makedirs(log_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise e

    logger = logging.getLogger(__name__)

    log_file = "{}/{}_{}.log".format(log_dir, log_file_basename, log_file_suffix)

    try:
        logging_config.dictConfig(_create_config(log_level, log_file))
        success_msg = "Logging configuration was loaded. Log messages can be found at {}.".format(log_file)
        logger.info(success_msg)
    except Exception as e:
        logger.error("Failed to load logging config!")
        raise e


def _create_config(log_level, log_file):
    return {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'simple': {
                'format': '%(asctime)s - %(module)s:%(lineno)d - %(levelname)s - %(message)s',
            },
            'detailed': {
                'format': '%(name)s:%(levelname)s %(module)s:%(lineno)d:  %(message)s',
            },
        },
        'handlers': {
            'console': {
                'level': log_level,
                'class': 'logging.StreamHandler',
                'formatter': 'simple',
                'stream': sys.stdout,
            },
            'file': {
                'level': log_level,
                'class': 'logging.FileHandler',
                'formatter': 'simple',
                'filename': log_file,
                'mode': 'w',
            },
        },
        'loggers': {
            '': {
                'handlers': ['console', 'file'],
                'level': log_level,
            },
        },
    }
