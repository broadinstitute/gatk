Python packages to be used in the GATK should be placed in this directory.  
Each package should be contained in its own subdirectory.  If the package 
can be installed as a standalone package, a corresponding `setup_<PACKAGE_NAME>.py` 
file may be placed in this directory.  However, during creation of the common 
GATK conda environment, all packages will be combined and pip-installed as a 
single package named ``gatkpythonpackages`` by `setup.py`.
