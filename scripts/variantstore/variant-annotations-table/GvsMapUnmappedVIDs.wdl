version 1.0

import "../wdl/GvsUtils.wdl" as Utils

workflow GvsMapUnmappedVIDs {
    input {
        File sites_only_vcf
        String project
        String dataset
        String vat_table_name
        String unmapped_vid_mapping_table_name
        String reference_name = "hg38"
    }

    call Utils.GetToolVersions

    call Utils.GetReference {
        input:
            reference_name = reference_name,
            basic_docker = GetToolVersions.basic_docker,
    }

    call MapUnmappedVIDs {
        input:
        project = project,
        dataset = dataset,
        vat_table_name = vat_table_name,
        sites_only_vcf = sites_only_vcf,
        reference = GetReference.reference,
        unmapped_vid_mapping_table_name = unmapped_vid_mapping_table_name,
        variants_docker = GetToolVersions.variants_docker,
    }

}

task MapUnmappedVIDs {
    input {
        String project
        String dataset
        String vat_table_name
        File sites_only_vcf
        Reference reference
        String unmapped_vid_mapping_table_name
        String variants_docker
    }
    parameter_meta {
        sites_only_vcf: {
            localization_optional: true
        }
    }
    command <<<

        # 1. Get the unmapped VIDs.
        bq --apilog=false query --use_legacy_sql=false --project_id=~{project} --format=csv --max_rows 1000000000 '

            select distinct vid from `~{dataset}.~{vat_table_name}` where vid not in (select vid from `~{dataset}.~{vat_table_name}`)

        ' sed 1d > unmapped_vids.csv

        # 2. Generate bcftools commands to query the sites-only VCF.
        python generate_bcftools_commands.py unmapped_vids.csv | sed 's/$/ >> to_search.vcf/' > bcftools_commands.sh

        # 3. Initialize a `to_search.vcf` with a VCF header, using the sites-only VCF as a donor.
        bcftools head sites_only.vcf > to_search.vcf

        # 4. Run the generated bcftools commands to extract the relevant variants from the sites only VCF and add the
        #    output to this `to_search.vcf` file.
        bash bcftools_commands.sh

        # 5. The generated bcftools commands will likely have overlapping ranges and generate duplicate and non-sequential entries
        #    in the `to_search.vcf` file. Fix that by sorting and deduplicating the VCF.
        bcftools sort to_search.vcf | bcftools norm -d none -o to_search.sort.dedup.vcf

        # 6. Make a normalized (left aligned) version of this VCF.
        bcftools norm -f ~{reference.reference_fasta} to_search.sort.dedup.vcf > to_search.sort.dedup.norm.vcf

        # 7. At this point the normalized VCF should contain all the "pseudo vids", but because the bcftools searches are written
        #    to search a range of positions, they will likely match some variants that do not correspond to "pseudo vids".
        #    Run the following to clean out non-"pseudo vid" entries from this file.
        python filter_vcf_by_vids.py pseudo_vids_file.tsv to_search.sort.dedup.norm.vcf > hits_only.vcf

        # 8. Now we can correlate the entries in this `hits_only.vcf` file back to the non-left aligned version that uses the same
        #    positions as GVS.
        python compare_vcfs.py to_search.sort.dedup.vcf hits_only.vcf > unmapped_vid_mappings.tsv

        # 9. Load into BigQuery
        bq load --project_id ~{project} --source_format=CSV --skip_leading_rows=1 --field_delimiter="\t" \
        ~{dataset}.~{unmapped_vid_mapping_table_name} unmapped_vid_mappings.tsv \
        vid:STRING,chr:STRING,input_location:INTEGER,input_position:INTEGER,input_ref:STRING,input_alt:STRING,left_aligned_location:INTEGER,left_aligned_position:INTEGER,left_aligned_ref:STRING,left_aligned_alt:STRING,info_field:STRING

    >>>
    runtime {
        docker: variants_docker
    }
}
