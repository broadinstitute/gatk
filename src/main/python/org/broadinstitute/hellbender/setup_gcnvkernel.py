from distutils.core import setup
import re
import sys

assert sys.version_info >= (3, 4), "gcnvkernel requires Python 3.4.x or later"

def get_version_string():
    version_file = "gcnvkernel/_version.py"
    version_str_line = open(version_file, "rt").read()
    version_regexp = r"^__version__ = ['\"]([^'\"]*)['\"]"
    re_out = re.search(version_regexp, version_str_line, re.M)
    if re_out is not None:
        return re_out.group(1)
    else:
        raise RuntimeError("Unable to find version string in %s." % (version_file,))

setup(
    name='gcnvkernel',
    version=get_version_string(),
    author='Mehrtash Babadi',
    author_email='mehrtash@broadinstitute.org',
    packages=['gcnvkernel',
              'gcnvkernel.inference',
              'gcnvkernel.models',
              'gcnvkernel.preprocess',
              'gcnvkernel.postprocess',
              'gcnvkernel.structs',
              'gcnvkernel.tasks',
              'gcnvkernel.utils',
              'gcnvkernel.io'],
    license='LICENSE.txt',
    description='GATK gCNV computational kernel',
    long_description=open('gcnvkernel/README.txt').read(),
    install_requires=[
        "theano == 0.9.0",
        "pymc3 == 3.1",
        "numpy >= 1.13.1",
        "scipy >= 0.19.1",
        "tqdm >= 4.15.0"
    ])
