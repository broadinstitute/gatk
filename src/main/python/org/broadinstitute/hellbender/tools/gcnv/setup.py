from distutils.core import setup
import re
VERSIONFILE="gcnvkernel/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(
    name='gcnvkernel',
    version=verstr,
    author='Mehrtash Babadi',
    author_email='mehrtash@broadinstitute.org',
    packages=['gcnvkernel',
              'gcnvkernel.inference',
              'gcnvkernel.models',
              'gcnvkernel.preprocess',
              'gcnvkernel.structs',
              'gcnvkernel.tasks',
              'gcnvkernel.utils',
              'gcnvkernel.io'],
    license='LICENSE.txt',
    description='GATK gCNV computational kernel',
    long_description=open('README.txt').read(),
    install_requires=[
        "theano == 0.9.0",
        "pymc3 == 3.1",
        "numpy >= 1.13.1",
        "scipy >= 0.19.1",
        "tqdm >= 4.15.0",
        "argparse"
    ],
)
