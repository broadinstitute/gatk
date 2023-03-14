# Creating a VDS for AoU

**NOTE:** This doc is very much a work in process and is serving as a way to keep track of Hail commands that were used for Delta after the initial VDS delivery.
This example walks through the procedure we followed for creating a VDS.

## Configuring the Terra notebook cluster

Please follow the instructions
in [AoU Delta VDS Cluster Configuration](cluster/AoU%20Delta%20VDS%20Cluster%20Configuration.md) to set up a Delta-scale
cluster.


## Now we made a VDS for Delta, and then needed to run the following to correct the GQ0s and then re-deliver it

```
pip install --force-reinstall  hail==0.2.109
python3

>>> import hail as hl
>>> hl.init(tmp_dir='gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/temp/')
>>> vds = hl.vds.read_vds('gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/vds/01-04-2023-correct-GT/aou_srwgs_short_variants_v7_without_ext_aian_prod.vds')
```

## Turn the GQ0s into no calls

```python
>>> rd = vds.reference_data
>>> rd = rd.filter_entries(rd.GQ > 0)
>>> rd = rd.filter_rows(hl.agg.count() > 0) # remove sites with no block starts
>>> rd = rd.annotate_globals(ref_block_max_length=1000) # set the max length
>>> vds2 = hl.vds.VariantDataset(reference_data=rd, variant_data=vds.variant_data)
```

## Put it back in the cloud

```python
>>> final_path = 'gs://prod-drc-broad/aou-wgs-delta/vds/2023-02-13-GQ0/aou_srwgs_short_variants_v7_without_ext_aian_prod_gq0.vds'
>>> vds2.write(final_path, overwrite=True)
```

## Validate it to ensure that it is ready to be shared

Copy the [VDS Validation python script](vds_validation.py) to the notebook environment.
Run it with the following arguments:

`--vds-path`: the GCS path to the newly-created VDS

`--temp-path`: a GCS path for temp files to live


An example run (on the Anvil 3k) is illustrated below

```
python3 validation_script.py --vds-path 'gs://fc-3ecd704c-4d2c-4360-bae1-093e214abce2/vds/vds_GQ0_max_ref_1k' --temp-path 'gs://fc-62891290-7fa3-434d-a39e-c64eeee4db8d/temp'
```