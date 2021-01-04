from distutils.core import setup
import re
import sys

assert sys.version_info >= (3, 4), "svgenotyper requires Python 3.4.x or later"

def get_version_string():
    version_file = "svgenotyper/_version.py"
    version_str_line = open(version_file, "rt").read()
    version_regexp = r"^__version__ = ['\"]([^'\"]*)['\"]"
    re_out = re.search(version_regexp, version_str_line, re.M)
    if re_out is not None:
        return re_out.group(1)
    else:
        raise RuntimeError("Unable to find version string in %s." % (version_file,))

setup(
    name='svgenotyper',
    version=get_version_string(),
    author='Mark Walker',
    author_email='markw@broadinstitute.org',
    packages=['svgenotyper'],
    license='LICENSE.txt',
    description='GATK structural variation genotyper',
    install_requires=[
        "pyro-ppl == 1.2.1"
    ])
