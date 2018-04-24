import re
from distutils.core import setup

def get_version_string():
    version_file = "vqsr_cnn/_version.py"
    version_str_line = open(version_file, "rt").read()
    version_regexp = r"^__version__ = ['\"]([^'\"]*)['\"]"
    re_out = re.search(version_regexp, version_str_line, re.M)
    if re_out is not None:
        return re_out.group(1)
    else:
        raise RuntimeError("Unable to find version string in %s." % (version_file,))

setup(name='vqsr_cnn',
      version=get_version_string(),
      description='Variant quality score recalibration with Convolutional Neural Networks',
      author='Sam Friedman',
      author_email='sam@broadinstitute.org',
      license='LICENSE.txt',
      packages=['vqsr_cnn'],
      install_requires=[
          "keras >= 2.0",
          "numpy >= 1.13.1",
          "scipy >= 0.19.1",
          "pysam >= 0.13",
          "scikit-learn >= 0.19.1",
          "matplotlib >= 2.1.2",
          "pyvcf >= 0.6.8",
          "biopython >= 1.70"
      ]
)
