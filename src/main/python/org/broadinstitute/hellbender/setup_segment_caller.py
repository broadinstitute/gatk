from distutils.core import setup
import re

def get_version_string():
    version_file = "segment_caller/_version.py"
    version_str_line = open(version_file, "rt").read()
    version_regexp = r"^__version__ = ['\"]([^'\"]*)['\"]"
    re_out = re.search(version_regexp, version_str_line, re.M)
    if re_out is not None:
        return re_out.group(1)
    else:
        raise RuntimeError("Unable to find version string in %s." % (version_file,))

setup(
    name='segment_caller',
    version=get_version_string(),
    author='Marton Kanasz-Nagy',
    author_email='mkanaszn@broadinstitute.org',
    packages=['segment_caller'],
    description='Calls copy number variation events on each segment, taking '
                'into account both the copy ratio and the allele fraction data.',
    long_description=open('gcnvkernel/README.txt').read(),
    install_requires=[
        "scikit-learn >= 0.18.1",
        "scipy >= 0.19.0",
        "matplotlib >= 2.1.0",
        "numpy >= 1.13.1"
    ])
