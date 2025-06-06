# Introduction

The following is a procedure for mapping the "pseudo vids" in the VAT that don't connect to `alt_allele` or `filter_info` 
entries in GVS. The `.prompt` files in this directory were used with Claude Sonnet 4 to generate the similarly named
sibling scripts. Please consult these files for more details on each script.

# Pseudo VID resolution procedure on a Terra notebook terminal

## Gather required data
  - `gcloud storage cp` the sites only VCF and index used as input to the VAT. The sites only is huge and takes a while to download, kick it off in another `screen` window.
 
  - `gcloud storage cp` the HG38 reference FASTA and index files as specified in `GvsUtils.GetReference`. These are also pretty big and good candidates for another `screen` window.
 
  - Download the "pseudo vid" data. For Echo this is in `aou-genomics-curation-prod.echo.pseudo_vids_only_in_vat`:

    ```shell
    bq --apilog=false query --max_rows 10000000 --project_id aou-genomics-curation-prod --use_legacy_sql=false --format=csv '
        SELECT vid FROM `aou-genomics-curation-prod.echo.pseudo_vids_only_in_vat`
    ' | sed 1d > pseudo_vids_file.tsv
    ```

## Select out the relevant parts of the sites only VCF
  - From manually examining a few of the smaller "pseudo vids" a pattern emerged of synonymous, normalized but non-left aligned variants in the sites-only VCF.
    With a few iterations, a script for generating bcftools commands to select out the relevant portions of the sites only VCF was created that
    identified non-left aligned correspondences for all the "pseudo vids" in the VAT.

    The `generate_bcftools_commands.py` script generates bcftools commands to select out relevant portions of the sites only VCF (see prompt for details):
 
    ```shell
    python generate_bcftools_commands.py pseudo_vids_file.tsv | sed 's/$/ >> to_search.vcf/' > bcftools_commands.sh
    ```

  - Initialize a "to_search.vcf" with a VCF header. Probably any of our hg38 VCFs will work, but the sites-only VCF should be handy:
    ```shell
    bcftools head sites_only.vcf > to_search.vcf
    ```
  - Run the generated bcftools commands to extract the relevant variants from the sites only VCF and add the to this "to_search.vcf" file:
    ```shell
    bash bcftools_commands.sh
    ```

  - The generated bcftools commands will likely have overlapping ranges and generate duplicate and non-sequential entries
    in the "to_search.vcf" file. Fix that by sorting and deduplicating the VCF:
    ```shell
    bcftools sort to_search.vcf | bcftools norm -d none -o to_search.sort.dedup.vcf
    ```
    
## Normalize and filter

  - Now make a normalized (left aligned) version of this VCF:
    ```shell
    bcftools norm -m- -f Homo_sapiens_assembly38.fasta to_search.sort.dedup.vcf > to_search.sort.dedup.norm.vcf
    ```

  - At this point the normalized VCF should contain all the "pseudo vids", but because the bcftools searches are not
    overly specific, they will also match some variants that do not correspond to "pseudo vids". Run the following
    command to clean out non-"pseudo vid" entries from this file:
    ```shell
    python filter_vcf_by_vids.py pseudo_vids_file.tsv to_search.sort.dedup.norm.vcf > hits_only.vcf
    ```
    
## Correlate left aligned to input positions, load into BigQuery
  - Now we can correlate the entries in this "hits_only.vcf" file back to the non-left aligned version that uses the same
    positions as GVS (see the prompt and code for details of how this works):
    ```shell
    python compare_vcfs.py to_search.sort.dedup.vcf hits_only.vcf > pseudo_vid_mappings.tsv
    ```
  - Load into BigQuery:
    ```shell
    bq load --project_id aou-genomics-curation-prod --source_format=CSV --skip_leading_rows=1 --field_delimiter="\t" \
        echo.pseudo_vid_mapping pseudo_vid_mappings.tsv \
        vid:STRING,chr:STRING,input_location:INTEGER,input_position:INTEGER,input_ref:STRING,input_alt:STRING,left_aligned_location:INTEGER,left_aligned_position:INTEGER,left_aligned_ref:STRING,left_aligned_alt:STRING,info_field:STRING
    ```