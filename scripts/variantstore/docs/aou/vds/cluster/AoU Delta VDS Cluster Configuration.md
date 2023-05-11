# Delta VDS Hail Cluster Configuration

To work with a Delta-sized VDS (roughly 250K samples) we use our PMI ops accounts to create a cluster in a Terra
customized Jupyter notebook with the following configuration:

![AoU Delta Cluster Config](AoU%20Delta%20VDS%20Cluster%20Configuration.png)

The startup script can be found at the following location. Permissions may need to be adjust to allow one's
Terra proxy group to read this object:

```
gs://gvs-internal/bigquery-jointcalling/hail/delta/patch_memory_2022_11_17.sh
```

The script has the contents:

```shell
#!/bin/sh
printf "\nspark.driver.memory=85g\n" >> /etc/spark/conf/spark-defaults.conf
```

**Note:** this script hardcodes a memory value for the driver node. If using a driver node with a different amount of
memory this value should be adjusted to avoid using too little or too much memory. Please ask the Hail team in Zulip if
uncertain what value to use.

## Installing the current Hail wheel

Once the cluster is up and running, open a notebook terminal. Within the terminal copy down the latest Hail wheel. **
Note:** you will need the proxy group of your PMI ops account to have at least the "Storage Object Viewer" role on the
bucket where the Hail wheel lives.

```
$ gsutil -m cp gs://gvs-internal-scratch/hail-wheels/2022-10-18/0.2.102-964bee061eb0/hail-0.2.102-py3-none-any.whl .
```

Install the Hail wheel and its dependencies. At the time of this writing the installation process produces a lot of
warnings and messages in red text which appear to be benign and can be ignored.

```
$ pip install --force-reinstall hail-0.2.102-py3-none-any.whl
```

## Start a Python REPL and import Hail

Finally start up a Python REPL and import Hail:

```python
>>> import hail as hl
```

You're now ready to work with a Delta-scale Hail VDS. 
