version 1.0

import "../wdl/GvsUtils.wdl" as Utils
import "../structs/Range.wdl" as Range

workflow GvsMapDroppedDuplicateVIDs {
    input {
        File sites_only_vcf
        File filtered_synonyms
        String project_id
        String dataset
        String participant_mapping_table_name
        String duplicate_mapping_table_name
        String reference_name = "hg38"
        Range? range_filter
    }

    call Utils.GetToolVersions

    call Utils.GetReference {
        input:
            reference_name = reference_name,
            basic_docker = GetToolVersions.basic_docker,
    }

    # Intentionally unused: this workflow exists solely to run this task as a side effect; there are no outputs to consume.
    #@ except: UnusedCall
    call MapDroppedDuplicateVIDs {
        input:
            project_id = project_id,
            dataset = dataset,
            participant_mapping_table_name = participant_mapping_table_name,
            sites_only_vcf = sites_only_vcf,
            filtered_synonyms = filtered_synonyms,
            reference = GetReference.reference,
            duplicate_mapping_table_name = duplicate_mapping_table_name,
            range_filter = range_filter,
            variants_docker = GetToolVersions.variants_docker,
    }
}

task MapDroppedDuplicateVIDs {
    input {
        String project_id
        String dataset
        String participant_mapping_table_name
        File sites_only_vcf
        File filtered_synonyms
        Reference reference
        String duplicate_mapping_table_name
        Range? range_filter
        String variants_docker
    }
    parameter_meta {
        sites_only_vcf: {
            localization_optional: true
        }
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # 1. auth before doing GCS-flavored bcftools commands
        set +o xtrace
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        set -o xtrace

        # 2. Subsequent steps require an index for the sites-only VCF. If an index is already present in the cloud just
        #    use that, otherwise download the sites-only VCF and index it.
        gsutil stat ~{sites_only_vcf}.tbi

        if [ $? -ne 0 ]
        then
            gcloud storage cp ~{sites_only_vcf} .
            sites_only="$(basename ~{sites_only_vcf})"
            bcftools index --tbi ${sites_only}
        else
            sites_only="~{sites_only_vcf}"
        fi

        # 3. Get the dropped duplicate VCFs. Transform the filtered synonyms file into an actual VCF by adding a header,
        #    kindly donated by the sites-only VCF.
        bcftools head ${sites_only} > filtered_synonyms.vcf
        cp filtered_synonyms.vcf head.vcf
        cat ~{filtered_synonyms} >> filtered_synonyms.vcf

        # 4. Get unique VIDs of the dropped duplicates.
        # Process the autosomal chromosomes first. The X and Y chromosomes complicate sorting; add in X and Y in a
        # separate step, where the chromosome field is not sorted numerically. The resulting file will have chromosomes
        # sorted in numerical order followed by X, followed by Y; secondary sort on numerical position; tertiary sort
        # on reference allele; quaternary sort on alt allele.
        bcftools query -f '%CHROM-%POS-%REF-%ALT\n' filtered_synonyms.vcf | cut -c 4- | grep -E -v '^[XY]' | sort -t - -k 1,1n -k 2,2n -k 3,3 -k 4,4 | uniq > dropped_duplicate_vids.tsv
        bcftools query -f '%CHROM-%POS-%REF-%ALT\n' filtered_synonyms.vcf | cut -c 4- | grep -E '^[XY]' | sort -t - -k 1,1 -k 2,2n -k 3,3 -k 4,4 | uniq >> dropped_duplicate_vids.tsv

        # 5. Generate bcftools commands to query the sites-only VCF.
        python /app/generate_bcftools_searches_for_variant_synonyms.py --sites-only ${sites_only} dropped_duplicate_vids.tsv | sed 's/$/ >> to_search.vcf/' > bcftools_commands.sh

        # 6. Initialize a `to_search.vcf` with a VCF header, using the sites-only VCF as a donor.
        bcftools head ${sites_only} > to_search.vcf

        # 7. Run the generated bcftools commands to extract the relevant variants from the sites only VCF and add the
        #    output to this `to_search.vcf` file.
        bash bcftools_commands.sh

        # 8. The generated bcftools commands will likely have overlapping ranges and generate duplicate and non-sequential entries
        #    in the `to_search.vcf` file. Fix that by sorting and deduplicating the VCF.
        bcftools sort to_search.vcf | bcftools norm -d none -o to_search.sort.dedup.vcf

        # 9. Make a normalized (left aligned) version of this VCF.
        bcftools norm -f ~{reference.reference_fasta} to_search.sort.dedup.vcf > to_search.sort.dedup.norm.vcf

        # 10. At this point the normalized VCF should contain all the duplicate vids, but because the bcftools searches are written
        #     to search a range of positions, they will likely also match some variants that do not correspond to duplicate vids.
        #     Run the following to clean out non-duplicate vid entries from this file.
        python /app/filter_vcf_by_vids.py dropped_duplicate_vids.tsv to_search.sort.dedup.norm.vcf > hits_only.vcf

        # 11. Now we can correlate the entries in this `hits_only.vcf` file back to the non-left aligned version that uses the same
        #     positions as GVS.
        python /app/map_input_alignments_to_left_alignments.py to_search.sort.dedup.vcf hits_only.vcf > dropped_duplicate_mappings.tsv

        # 12. Load the duplicate vid mapping into BigQuery
        bq load --project_id ~{project_id} --source_format=CSV --skip_leading_rows=1 --field_delimiter="\t" \
            ~{dataset}.~{duplicate_mapping_table_name} dropped_duplicate_mappings.tsv \
            vid:STRING,chr:STRING,input_location:INTEGER,input_position:INTEGER,input_ref:STRING,input_alt:STRING,left_aligned_location:INTEGER,left_aligned_position:INTEGER,left_aligned_ref:STRING,left_aligned_alt:STRING,info_field:STRING

        # 13. Remove existing mappings for these duplicates from the mapping table
        #
        # bq query --max_rows check: ok delete
        #
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} --format=csv '

        DELETE `~{dataset}.~{participant_mapping_table_name}` WHERE vid IN (
            SELECT dup.vid FROM `~{dataset}.~{duplicate_mapping_table_name}` dup
        ~{if defined(range_filter) then "WHERE dup.left_aligned_location >= ~{select_first([range_filter]).startLocation} AND dup.left_aligned_location < ~{select_first([range_filter]).endLocation}" else ""}
        )

        '

        # 14. Write all the person mappings for these duplicate VIDs
        #
        # bq query --max_rows check: ok insert
        #
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} --format=csv '

        INSERT into `~{dataset}.~{participant_mapping_table_name}` (vid, person_ids)
        SELECT
          dup.vid AS vid,
          ARRAY_AGG(SAFE_CAST(si.sample_name AS INT64) IGNORE NULLS) AS person_ids
        FROM
          `~{dataset}.~{duplicate_mapping_table_name}` dup
        JOIN
          `~{dataset}.alt_allele` aa
        ON
          dup.input_location = aa.location
          AND dup.input_ref = aa.ref
          AND dup.input_alt = aa.allele
        JOIN
          `~{dataset}.sample_info` si
        ON
          aa.sample_id = si.sample_id
          -- `alt_allele` is append-only and is only ever extended for sample ids greater than the maximum
          -- already present, so it still contains rows for samples that were withdrawn after those rows were
          -- written. Those samples are excluded from the VDS, so they must be excluded here too. Controls are
          -- likewise absent from the VDS; the SAFE_CAST above happens to drop non-numeric control sample names,
          -- but do not rely on that.
          AND si.withdrawn IS NULL
          AND si.is_control = false
          ~{if defined(range_filter) then "WHERE aa.location >= ~{select_first([range_filter]).startLocation} AND aa.location < ~{select_first([range_filter]).endLocation}" else ""}
        GROUP BY
          vid

        '
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Boolean done = true
        Array[File] shell_scripts = glob("*.sh")
        Array[File] vcfs = glob("*.vcf")
        Array[File] tsvs = glob("*.tsv")
    }
}

