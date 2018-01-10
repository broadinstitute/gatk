# Taken from http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
import os


def count_lines(filename):
    """
    returns integer of number of lines
    :param filename:  name of file.
    :type filename: str
    :return: number of lines in the file
    """
    if not os.path.exists(filename):
        raise IOError("Please make sure that " + filename + " exists.  Could not be found or read.")
    lines = 0
    with open(filename) as f:
        buf_size = 1024 * 1024
        read_f = f.read # loop optimization

        buf = read_f(buf_size)
        while buf:
            lines += buf.count('\n')
            buf = read_f(buf_size)

    return lines

