# Callset cleanup

During the normal course of creating an AoU callset several large and expensive artifacts are created:

* BigQuery dataset
* Avro files (Delta onward)
* Hail VariantDataset (VDS) (Delta onward)

The current Variants AoU callset policy is effectively to retain all of these artifacts for all callset versions
forever. As there can be very significant storage costs involved, the Variants team would like to make artifacts costs
more clear so AoU can make conscious choices about which artifacts they actually want to keep.

## Establish an AoU callset signoff protocol

### Existing internal signoff protocol

The Variants team currently has the following internal VDS signoff protocol:

* Generate a VDS for the candidate callset
* Generate callset statistics for the candidate callset
* Forward both items above to Lee for QA / approval

Successful completion of this internal signoff protocol gates
the [delivery of the VDS to AoU](../vds/delivery/Delivering%20a%20VDS.md).

### Future AoU signoff protocol

Analogous to the existing internal signoff protocol, the Variants team would like to have an external AoU callset
signoff protocol. The completion of this protocol could gate removal of intermediate data such as Avro files and
pre-delivery copies of the VDS.

There are tradeoffs in cost and time between retaining intermediate callset data versus regenerating it if required,
details of which are discussed below. Ultimately AoU is the party paying for these cloud resources so AoU should
decide what they would like the cleanup policies to be.

### Avro file storage versus regeneration

The Avro files generated from the Delta callset onward are very large, significantly larger than the final Hail VDS.
For the ~250K sample Delta callset the Avro files consume nearly 80 TiB of GCS storage while the delivered VDS is
"only" about 26 TiB.

If the Avro files needed to be regenerated there would be a significant cost to re-scan the BigQuery tables that
were used to make them. On the other hand, continuing to store the existing and currently unused Avro files incurs
significant GCS storage fees. Per
the [cost spreadsheet](https://docs.google.com/spreadsheets/d/1fcmEVWvjsx4XFLT9ZUsruUznnlB94xKgDIIyCGu6ryQ/edit#gid=0),
for the Delta callset, recreating the Avro files would take about 12 hours and cost approximately $3000 in BigQuery scan
charges. Alternatively, continued storage of Delta Avro files costs approximately $1500 / month.

See "AoU signoff protocol" in [Tickets to be created](#tickets-to-be-created)below.

### Hail VariantDataset (VDS) storage versus regeneration

The Hail VDS generated for the Delta callset consumes about 26 TiB of space in GCS at a cost of approximately $500 /
month. Recreating the VDS from Avro files would take around 10 hours at about $100 / hour in cluster time for a total of
about $1000. Note that re-creating the VDS requires Avro files; if we have not retained the Avro files per the step
above, we would need to regenerate those as well which would add significantly to the cost.

See "AoU signoff protocol" in [Tickets to be created](#tickets-to-be-created)below.

## "Inflection points"

With each new AoU callset the Variants team has pushed the limits of scale beyond any previous callsets
we’ve created. Many times this requires us to plan for new approaches where we know we’ve hit scale limitations, but
other times the team only discovers that we’ve hit a scale limitation once the callset is underway.

While a callset is underway the team might receive a request from Lee and/or AoU to make changes to the format,
contents, or delivery of callset data. These requests may invalidate some of the data that has been produced up to that
point in the callset.

And sometimes we might just mess something up.

These "inflection points" in callset creation should motivate the Variants team to think about what data generated
up to that point should be deleted. This should not only help make clear the work that does belong to the current
callset, but would also prevent spending money on the storage of data that is not part of the callset.

Some of this data can consume significant space: for Delta we found ~300 TiB of junk costing several thousand dollars
per month. Cleanup usually includes separately deleting the BQ dataset and workspace as the deletion of one does not
cascade to the other.

For the AoU Delta callset, the Variants team delivered the callset in Hail VDS format for the first time. The team’s
lack of familiarity with Hail unsurprisingly caused us to make several wrong turns trying to make this callset:

1. Our first attempt to extract Avro files produced unusable results since an incorrect filter name was supplied and the
   code at the time did not validate filter names. We noticed the mistake fairly soon after completing the export, but
   unfortunately it took a while for us to realize that the unusable Avro output was consuming ~80 TiB of GCS storage.
2. Some of our intermediate callsets used a GQ 40 drop state which turned out to be unworkable for Hail. Once we
   realized this and committed to a "NONE" drop state we should have promptly cleaned up unusable BQ datasets and
   associated workspaces.
3. The Delta v1 workspace + BQ dataset was similarly abandoned due to improperly reblocked VCFs and should have been
   deleted promptly.
4. We unexpectedly found ourselves running the "prepare" step in order to generate callset statistics as we had not
   realized during callset planning that callset statistics generation depended on the tables created and populated in
   this step. While we only needed the variant data to generate callset statistics, at the time the "prepare" workflow
   did not have the ability to bypass generation of reference and sample data tables. The "prepare" reference data in
   particular consumed ~30 TiB and should have been deleted right away.
5. The Hail GVS import script specifies a temporary "directory” which is not automatically cleaned up if the import
   fails. This was useful for our team when the Hail import was failing after creating intermediate VDSes since we were
   able to pick up where we left off on subsequent attempts and save a lot of compute. However once the import finally
   succeeded we did not immediately realize that the temporary directory is still not cleaned up automatically, and
   furthermore that this directory’s contents are several times larger than the actual output VDS.

## Historical cleanup

### BigQuery datasets

Several large datasets exist within the `aou-genomics-curation-prod` project, most of
which would now be considered historical (methodology [here](#code-for-querying-datasets-and-their-sizes)):

| Name | Size (TiB) | Notes              |
|------|------------|--------------------|
|alpha1_1000|6.17| ?                  |
|alpha2_prelim2|205.41| ?                  |
|aou_wgs|469.92| Beta + Charlie?    |
|aou_wgs_10k|53.18| ?                  |
|aou_wgs_fullref_v2|691.2| Delta              |
|batch_effects_test_v1|15.25| ?                  |
|beta2_99k|385.84| Beta intermediate? |
|beta_release|1004.57| Beta?              |
|cdr_metadata|0.0| ?                  |
|nosa_wdl_test_v1|0.06| ?                  |
|rc_add_AD_1000|3.84| ?                  |
|temp_tables|0.0| ?                  |

# Appendix

## Tickets to be created

* AoU signoff protocol
    * Talk to Lee to identify people on the AoU side who could sign off on a delivered callset
    * Once the AoU decisionmakers have been identified, explain the intermediate data we're currently keeping and the
      cost tradeoffs of keeping versus regenerating data. In particular note that for Delta we are currently keeping two
      copies of the VDS (though not two sets of Avro files).
* Identifying existing BigQuery datasets
    * Document the purpose of the existing AoU datasets
    * Classify these datasets by: "Should be deleted", "Should not be deleted", and "Unsure".
    * Make any tickets required to attain certainty over whether datasets should be retained or not.
        * Can we copy filters from previous callsets into the awkwardly named Delta callset?
* Delete any "Should be deleted" datasets identified above.

## Code for querying datasets and their sizes

With PMI ops auth:

```shell
i=0
echo '[' > out.json
for dataset in $(bq ls --project_id aou-genomics-curation-prod --format prettyjson | jq -r '..|.datasetId? // empty')
do
    if [[ $i -ne 0 ]]; then
        echo ',' >> out.json
    fi
    i=$((i+1))

    bq query --nouse_legacy_sql --project_id=aou-genomics-curation-prod --format=prettyjson "
    
    SELECT
        '${dataset}' as dataset,
        ROUND(SUM(size_bytes) / POWER(2, 40), 2) as tebibytes
    FROM
      \`aou-genomics-curation-prod.${dataset}.__TABLES__\`

    " >> out.json
done
echo ']' >> out.json

cat out.json | jq -r '.[] | .[0] | [.dataset, .tebibytes] | @tsv' | sed "s/\t/|/" | sed -E 's/(.*)/|\1|/'
```
