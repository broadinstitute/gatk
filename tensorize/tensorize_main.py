import logging

from tensorize import tensorize
from tensorize.defines import REQUIREMENTS_FILE, SETUP_FILE


if __name__ == "__main__":
    logging.getLogger().setLevel(logging.DEBUG)
    tensorize.run(['--requirements_file={}'.format(REQUIREMENTS_FILE), '--setup_file={}'.format(SETUP_FILE)])
