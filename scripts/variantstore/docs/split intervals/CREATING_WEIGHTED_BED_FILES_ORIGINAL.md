Original issue link: https://github.com/broadinstitute/gatk/issues/7622


UPDATE: @lbergelson and I started creating this as an extension to SplitIntervals, but it quickly because very complex to fit it into that framework/abstraction so we decided to create a specialized tool for GVS

It would be valuable for SplitIntervals to be able to split intervals based not on number of genomic bases, but by using a set of weights.

Ideally this new mode would read in a BED file containing the weights in the score field and attempt to produce a series of intervals that have equal total weights.

Note: ` --dont-mix-contigs` should still continue to work

** Why? **
In the Genomic Variant Store, we have found that scattering work by "# of genomic bases" does not lead to even runtimes for the shards.

![image](https://user-images.githubusercontent.com/1423491/147964102-d2c83dea-8486-4699-9e7a-eef3dc759732.png)

Instead we have found that an excellent proxy for runtime is the number of variants contained in a given interval:

![image](https://user-images.githubusercontent.com/1423491/147964259-333f4058-a701-4b31-b410-242518d1b3b2.png)

And furthermore, that this generalizes even when we use a subset of a different dataset

![image](https://user-images.githubusercontent.com/1423491/147964392-17d045e9-e2f2-467b-8eae-77bd63291902.png)


The weights can be generated from a BQ data set like:

```
CREATE OR REPLACE TABLE `example.mydataset.vet_weight_100k` AS
SELECT CAST(TRUNC(location / 100000) * 100000 AS INT64) bin, count(*) entries
FROM `example.mydataset.vet_001`
GROUP BY 1 ORDER BY 1;
```

and then converted to bed with this python:

```
import pandas as pd

location_offset = 1000000000000
binsize_kb = 100
infile = f"40K_vet_weight_{binsize_kb}k.csv"

w = pd.read_csv(infile)
w['contig'] = "chr" + (w['bin'].astype(int) / location_offset).astype(int).astype(str).str.replace("23","X").replace("24","Y")
w['start_position'] = w['bin'].astype(int) - (w['bin'] / location_offset).astype(int) * location_offset
w['end_position'] = w['start_position'] + binsize_kb*1000
w['name'] = "."
o = w[['contig', 'start_position','end_position', 'name', 'entries' ]]

o.to_csv(f"gvs_vet_weights_{binsize_kb}kb.bed",sep='\t',index=False,header=False)
```
UPDATE: Latest version of the code is in the python script convert_weight_file_CSV_to_bed.py

After this step, run the generated bed file through close_bed_file_gaps.py.  It will create 0 weight blocks in all of the areas where there were no samples, which is necessary for the new Java code.


I've generated 3 BED files for 100kb, 10kb and 1kb and placed them into

```
929.42 KiB  2022-01-03T18:56:48Z  gs://broad-dsp-spec-ops/gvs/weights/gvs_vet_weights_100kb.bed
  8.78 MiB  2022-01-03T18:56:50Z  gs://broad-dsp-spec-ops/gvs/weights/gvs_vet_weights_10kb.bed
 84.65 MiB  2022-01-03T18:56:56Z  gs://broad-dsp-spec-ops/gvs/weights/gvs_vet_weights_1kb.bed
```

FWIW -- the intervals we're trying to divide are `gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list`

UPDATE: After the changes to how WeightedSplitIntervals operates, the weighted bed file in the new format is located here: `gs://gvs_quickstart_storage/weights/gvs_full_vet_weights_1kb_padded.bed`


