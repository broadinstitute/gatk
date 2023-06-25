version 1.0

import "GvsUtils.wdl" as Utils
import "../variant_annotations_table/GvsCreateVATFilesFromBigQuery.wdl" as GvsCreateVATFilesFromBigQuery


workflow GvsCreateVATfromVDS {
    input {
        File input_sites_only_vcf
        File ancestry_file

        String project_id
        String dataset_name
        String filter_set_name
        String? vat_version

        Int effective_scatter_count = 10

        String output_path
        Int? split_intervals_disk_size_override
        Int? split_intervals_mem_override
        Int? merge_vcfs_disk_size_override
    }

    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

    File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File reference_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    File reference_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

    call MakeSubpopulationFilesAndReadSchemaFiles {
        input:
            input_ancestry_file = ancestry_file
    }

    call Utils.IndexVcf {
        input:
            input_vcf = input_sites_only_vcf
    }

    call Utils.SplitIntervals {
        input:
            intervals = interval_list,
            ref_fasta = reference,
            ref_fai = reference_index,
            ref_dict = reference_dict,
            scatter_count = effective_scatter_count,
            output_gcs_dir = output_path + "intervals",
            split_intervals_disk_size_override = split_intervals_disk_size_override,
            split_intervals_mem_override = split_intervals_mem_override,
    }

    String sites_only_vcf_basename = basename(basename(input_sites_only_vcf, ".gz"), ".vcf")

    scatter(i in range(length(SplitIntervals.interval_files))) {
        String interval_file_basename = basename(SplitIntervals.interval_files[i], ".interval_list")
        String vcf_filename = interval_file_basename + "." + sites_only_vcf_basename

        call Utils.SelectVariants {
            input:
                input_vcf = input_sites_only_vcf,
                input_vcf_index = IndexVcf.output_vcf_index,
                interval_list = SplitIntervals.interval_files[i],
                output_basename = vcf_filename
        }

        call RemoveDuplicatesFromSitesOnlyVCF {
            input:
                sites_only_vcf = SelectVariants.output_vcf,
                ref = reference
        }

        call StripCustomAnnotationsFromSitesOnlyVCF {
            input:
                input_vcf = RemoveDuplicatesFromSitesOnlyVCF.output_vcf,
                custom_annotations_header = MakeSubpopulationFilesAndReadSchemaFiles.custom_annotations_template_file,
                output_vcf_name = "${vcf_filename}.unannotated.sites_only.vcf",
                output_custom_annotations_filename = "${vcf_filename}.custom_annotations.tsv"
        }

        ## Use Nirvana to annotate the sites-only VCF and include the AC/AN/AF calculations as custom annotations
        call AnnotateVCF {
            input:
                input_vcf = StripCustomAnnotationsFromSitesOnlyVCF.output_vcf,
                output_annotated_file_name = "${vcf_filename}_annotated",
                custom_annotations_file = StripCustomAnnotationsFromSitesOnlyVCF.output_custom_annotations_file,
        }

        call PrepVtAnnotationJson {
            input:
                positions_annotation_json = AnnotateVCF.positions_annotation_json,
                output_file_suffix = "${vcf_filename}.json.gz",
                output_path = output_path,
        }

        call PrepGenesAnnotationJson {
            input:
                genes_annotation_json = AnnotateVCF.genes_annotation_json,
                output_file_suffix = "${vcf_filename}.json.gz",
                output_path = output_path,
        }

    }

    call Utils.MergeTsvs {
        input:
            input_files = RemoveDuplicatesFromSitesOnlyVCF.track_dropped,
            output_file_name = "${sites_only_vcf_basename}.dropped.tsv"
    }

    call BigQueryLoadJson {
        input:
            nirvana_schema = MakeSubpopulationFilesAndReadSchemaFiles.vat_schema_json_file,
            vt_schema = MakeSubpopulationFilesAndReadSchemaFiles.variant_transcript_schema_json_file,
            genes_schema = MakeSubpopulationFilesAndReadSchemaFiles.genes_schema_json_file,
            project_id = project_id,
            dataset_name = dataset_name,
            output_path = output_path,
            filter_set_name = filter_set_name,
            vat_version = vat_version,
            prep_vt_json_done = PrepVtAnnotationJson.done,
            prep_genes_json_done = PrepGenesAnnotationJson.done
    }

    call GvsCreateVATFilesFromBigQuery.GvsCreateVATFilesFromBigQuery {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            vat_table_name = BigQueryLoadJson.vat_table_name,
            output_path = output_path,
            merge_vcfs_disk_size_override = merge_vcfs_disk_size_override,
            precondition_met = BigQueryLoadJson.done
   }

    output {
        File final_tsv_file = GvsCreateVATFilesFromBigQuery.final_tsv_file
        File dropped_sites_file = MergeTsvs.output_file
    }
}

################################################################################

task MakeSubpopulationFilesAndReadSchemaFiles {
    input {
        File input_ancestry_file

        String schema_filepath = "/data/variant_annotation_table/schema/"
        String vat_schema_json_filename = "vat_schema.json"
        String variant_transcript_schema_json_filename = "vt_schema.json"
        String genes_schema_json_filename = "genes_schema.json"
    }
    String output_ancestry_filename =  "ancestry_mapping.tsv"
    String custom_annotations_template_filename =  "custom_annotations_template.tsv"

    command <<<
        set -e

        cp ~{schema_filepath}* .

        ## the ancestry file is processed down to a simple mapping from sample to subpopulation
        python3 /app/extract_subpop.py \
            --input_path ~{input_ancestry_file} \
            --output_path ~{output_ancestry_filename} \
            --custom_annotations_template_path ~{custom_annotations_template_filename}
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-23-alpine"
        memory: "1 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File vat_schema_json_file = vat_schema_json_filename
        File variant_transcript_schema_json_file = variant_transcript_schema_json_filename
        File genes_schema_json_file = genes_schema_json_filename

        File ancestry_mapping_list = output_ancestry_filename
        File custom_annotations_template_file = custom_annotations_template_filename
    }
}


task StripCustomAnnotationsFromSitesOnlyVCF {
    input {
        File input_vcf
        File custom_annotations_header
        String output_vcf_name
        String output_custom_annotations_filename
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    Int disk_size = ceil((size(input_vcf, "GB") + size(custom_annotations_header, "GB")) * 4) + 100

    command <<<
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        python3 /app/strip_custom_annotations_from_sites_only_vcf.py \
        --input_vcf ~{input_vcf} \
        --input_custom_annotations_tsv ~{custom_annotations_header} \
        --output_vcf ~{output_vcf_name} \
        --output_custom_annotations_tsv ~{output_custom_annotations_filename}

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-23-alpine"
        memory: "7 GiB"
        cpu: "2"
        preemptible: 3
        disks: "local-disk " + disk_size + " HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File output_vcf = output_vcf_name
        File output_custom_annotations_file = output_custom_annotations_filename
        File monitoring_log = "monitoring.log"
    }
}


task RemoveDuplicatesFromSitesOnlyVCF {
    input {
        File sites_only_vcf
        File ref
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    Int disk_size = ceil(size(sites_only_vcf, "GB") * 5) + 100

    # separate multi-allelic sites into their own lines, remove deletions and filtered sites and make a sites only vcf
    # while extracting and calculating the an/ac/af & sc by subpopulation into a tsv
    command <<<
        set -e

        bash ~{monitoring_script} > monitoring.log &

        # custom function to prepend the current datetime to an echo statement
        echo_date () { echo "`date "+%Y/%m/%d %H:%M:%S"` $1"; }

        echo_date "VAT: Convert input to BCF format"
        bcftools convert --threads 4 -O b -o sites_only.bcf ~{sites_only_vcf}

        echo_date "VAT: Calculating number of sites with Ns"

        ## track the dropped variants with N's in the reference (Since Nirvana cant handle N as a base, drop them for now)
        bcftools view --threads 4 -i 'REF~"N"' -O u sites_only.bcf | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > track_dropped.tsv

        echo_date "VAT: filter out sites with N's in the reference AND sites with AC=0"
        ## NOTE: Sites that were filtered out because of AC=0 are not recorded in the 'track_dropped.tsv' file, but can be
        ##       determined by examining the sites-only VCF provided to this WDL.
        bcftools view --threads 4 -e 'REF~"N" || AC=0' -O b sites_only.bcf -o filtered_sites_only.bcf
        rm sites_only.bcf

        echo_date "VAT: normalize, left align and split multi allelic sites to new lines, remove duplicate lines"
        ## note that normalization may create sites with more than 50 alt alleles
        bcftools norm --threads 4 -m- --check-ref w -f ~{ref} filtered_sites_only.bcf -O b -o normalized.bcf
        rm filtered_sites_only.bcf

        echo_date "VAT: detecting and removing duplicate rows from sites-only VCF"

        ## During normalization, sometimes duplicate variants appear but with different calculations. This seems to be a bug in bcftools. For now we are dropping all duplicate variants
        ## to locate the duplicates, we first make a file of just the first 5 columns
        bcftools query normalized.bcf -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' | sort | uniq -d > duplicates.tsv

        echo_date "VAT: done with duplicate detection"
        wc -l duplicates.tsv
        echo_date "VAT: Duplicates may have been found"

        # If there ARE dupes to remove
        if [ -s duplicates.tsv ]; then
            ## remove those rows (that match up to the first 5 cols)
            echo_date "VAT: Removing those rows"
            bcftools view --threads 4 normalized.bcf | grep -v -wFf duplicates.tsv > deduplicated.vcf
        else
            # There are no duplicates to remove
            echo_date "VAT: No duplicates found"
            bcftools view --threads 4 normalized.bcf -o deduplicated.vcf
        fi
        rm normalized.bcf

        ## add duplicates to the file that's tracking dropped variants
        cat duplicates.tsv >> track_dropped.tsv
        rm duplicates.tsv ## clean up unneeded file

        echo_date "VAT: finished"
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-23-alpine"
        maxRetries: 3
        memory: "16 GB"
        preemptible: 3
        cpu: "8"
        disks: "local-disk " + disk_size + " HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File track_dropped = "track_dropped.tsv"
        File output_vcf = "deduplicated.vcf"
        File monitoring_log = "monitoring.log"
    }
}


task AnnotateVCF {
    input {
        File input_vcf
        String output_annotated_file_name
        File custom_annotations_file

        # Mentioning this path in the inputs section of the task combined with checking the 'Use reference disks' option
        # in Terra UI tells Cromwell to arrange for the Nirvana reference disk to be attached to this VM.
        File summon_reference_disk =
            "gs://broad-public-datasets/gvs/vat-annotations/Nirvana/3.18.1/SupplementaryAnnotation/GRCh38/MITOMAP_20200819.nsa.idx"
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String annotation_json_name = output_annotated_file_name + ".json.gz"
    String gene_annotation_json_name = output_annotated_file_name + ".genes.json.gz"
    String positions_annotation_json_name = output_annotated_file_name + ".positions.json.gz"
    String nirvana_location = "/Nirvana/Nirvana.dll"
    String custom_creation_location = "/Nirvana/SAUtils.dll"
    String jasix_location = "/Nirvana/Jasix.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"

    command <<<
        bash ~{monitoring_script} > monitoring.log &

        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # There's an issue with how the projects/broad-dsde-cromwell-dev/global/images/nirvana-3-18-1-references-2023-01-03
        # disk image was built: while all the reference files do exist on the image they are not at the expected
        # locations. The following code works around this issue and should continue to work even after a corrected
        # version of the Nirvana reference image is deployed into Terra.

        # Find where the reference disk should have been mounted on this VM.  Note this is referred to as a "candidate
        # mount point" because we do not actually confirm this is a reference disk until the following code block.
        CANDIDATE_MOUNT_POINT=$(lsblk | sed -E -n 's!.*(/mnt/[a-f0-9]+).*!\1!p')
        if [[ -z ${CANDIDATE_MOUNT_POINT} ]]; then
            >&2 echo "Could not find a mounted volume that looks like a reference disk, exiting."
            exit 1
        fi

        # Find one particular reference under the mount path. Note this is not the same reference as was specified in the
        # `inputs` section, so this would only be present if the volume we're looking at is in fact a reference disk.
        REFERENCE_FILE="Homo_sapiens.GRCh38.Nirvana.dat"
        REFERENCE_PATH=$(find ${CANDIDATE_MOUNT_POINT} -name "${REFERENCE_FILE}")
        if [[ -z ${REFERENCE_PATH} ]]; then
            >&2 echo "Could not find reference file '${REFERENCE_FILE}' under candidate reference disk mount point '${CANDIDATE_MOUNT_POINT}', exiting."
            exit 1
        fi

        # Take the parent of the parent directory of this file as root of the locally mounted references:
        DATA_SOURCES_FOLDER="$(dirname $(dirname ${REFERENCE_PATH}))"

        # =======================================
        echo "Creating custom annotations"
        mkdir customannotations_dir
        CUSTOM_ANNOTATIONS_FOLDER="$PWD/customannotations_dir"

        # Add AC/AN/AF as custom annotations
        ## use --skip-ref once you are on a version of nirvana later than 3.14 (once they have created a docker image for it)
        dotnet ~{custom_creation_location} customvar \
            -r $DATA_SOURCES_FOLDER~{path_reference} \
            -i ~{custom_annotations_file} \
            -o $CUSTOM_ANNOTATIONS_FOLDER

        # =======================================
        # Create Nirvana annotations:

        dotnet ~{nirvana_location} \
            -i ~{input_vcf} \
            -c $DATA_SOURCES_FOLDER~{path} \
            --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
            --sd $CUSTOM_ANNOTATIONS_FOLDER \
            -r $DATA_SOURCES_FOLDER~{path_reference} \
            -o ~{output_annotated_file_name}

        # https://illumina.github.io/NirvanaDocumentation/introduction/parsing-json#jasix
        # Parse out the Genes section into a separate annotated json
        dotnet  ~{jasix_location} \
            --in ~{annotation_json_name} \
            --section genes \
            --out ~{gene_annotation_json_name}

        # Parse out the Positions section into a separate annotated json
        dotnet  ~{jasix_location} \
        --in ~{annotation_json_name} \
        --section positions \
        --out ~{positions_annotation_json_name}

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:nirvana_2022_10_19"
        memory: "64 GB"
        cpu: "4"
        preemptible: 3
        disks: "local-disk 2000 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File genes_annotation_json = "~{gene_annotation_json_name}"
        File positions_annotation_json = "~{positions_annotation_json_name}"
        File monitoring_log = "monitoring.log"
    }
}

task PrepVtAnnotationJson {
    input {
        File positions_annotation_json
        String output_file_suffix
        String output_path
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String output_vt_json = "vat_vt_bq_load" + output_file_suffix
    String output_vt_gcp_path = output_path + 'vt/'

    command <<<
        set -o errexit -o nounset -o pipefail -o xtrace

        # Kick off the monitoring script
        bash ~{monitoring_script} > monitoring.log &

        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '

        ## the annotation jsons are split into the specific VAT schema
        python3 /app/create_vt_bqloadjson_from_annotations.py \
            --annotated_json ~{positions_annotation_json} \
            --output_vt_json ~{output_vt_json}

        gsutil cp ~{output_vt_json} '~{output_vt_gcp_path}'

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-23-alpine"
        memory: "7 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 500 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File vat_vt_json="~{output_vt_json}"
        Boolean done = true
        File monitoring_log = "monitoring.log"
    }
}

task PrepGenesAnnotationJson {
    input {
        File genes_annotation_json
        String output_file_suffix
        String output_path
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String output_genes_json = "vat_genes_bq_load" + output_file_suffix
    String output_genes_gcp_path = output_path + 'genes/'

    command <<<
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '

        ## the annotation jsons are split into the specific VAT schema
        python3 /app/create_genes_bqloadjson_from_annotations.py \
            --annotated_json ~{genes_annotation_json} \
            --output_genes_json ~{output_genes_json}

        gsutil cp ~{output_genes_json} '~{output_genes_gcp_path}'

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-23-alpine"
        memory: "7 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 500 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File vat_genes_json="~{output_genes_json}"
        Boolean done = true
        File monitoring_log = "monitoring.log"
    }
}


task BigQueryLoadJson {
    meta {
        # since the WDL will not see the updated data (its getting put in a gcp bucket)
        volatile: true
    }

    input {
        String filter_set_name
        String? vat_version
        File? nirvana_schema
        File? vt_schema
        File? genes_schema
        String project_id
        String dataset_name
        String output_path
        Array[Boolean] prep_vt_json_done
        Array[Boolean] prep_genes_json_done
    }

    # If the vat version is undefined or v1 then the vat tables would be named like filter_vat, otherwise filter_vat_v2.
    String effective_vat_version = if (defined(vat_version) && select_first([vat_version]) != "v1") then "_" + select_first([vat_version]) else ""

    # There are two pre-vat tables. A variant table and a genes table. They are joined together for the vat table
    String vat_table = filter_set_name + "_vat" + effective_vat_version
    String variant_transcript_table = filter_set_name + "_vat_vt" + effective_vat_version
    String genes_table = filter_set_name + "_vat_genes" + effective_vat_version

    String vt_path = output_path + 'vt/*'
    String genes_path = output_path + 'genes/*'

    command <<<
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        DATE=86400 ## 24 hours in seconds

        set +e
        bq --apilog=false show --project_id ~{project_id} ~{dataset_name}.~{variant_transcript_table} > /dev/null
        BQ_SHOW_RC=$?
        set -e

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
            bq --apilog=false --location=US mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{variant_transcript_table} ~{vt_schema}
        fi

        echo "Loading data into a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
        echo ~{vt_path}
        echo ~{genes_path}
        bq --apilog=false --location=US load --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON ~{dataset_name}.~{variant_transcript_table} ~{vt_path}

        set +e
        bq --apilog=false show --project_id ~{project_id} ~{dataset_name}.~{genes_table} > /dev/null
        BQ_SHOW_RC=$?
        set -e

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating a pre-vat table ~{dataset_name}.~{genes_table}"
            bq --apilog=false --location=US mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{genes_table} ~{genes_schema}
        fi

        echo "Loading data into a pre-vat table ~{dataset_name}.~{genes_table}"
        bq --apilog=false --location=US load  --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON  ~{dataset_name}.~{genes_table} ~{genes_path}

        set +e
        bq --apilog=false show --project_id ~{project_id} ~{dataset_name}.~{vat_table} > /dev/null
        BQ_SHOW_RC=$?
        set -e

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating the vat table ~{dataset_name}.~{vat_table}"
            bq --apilog=false --location=US mk --project_id=~{project_id} ~{dataset_name}.~{vat_table} ~{nirvana_schema}
        else
            bq --apilog=false rm -t -f --project_id=~{project_id} ~{dataset_name}.~{vat_table}
            bq --apilog=false --location=US mk --project_id=~{project_id} ~{dataset_name}.~{vat_table} ~{nirvana_schema}
        fi
        echo "And putting data into it"

        # Now we run a giant query in BQ to get this all in the right table and join the genes properly
        # Note the genes table join includes the group by to avoid the duplicates that get created from genes that span shards
        # Commented out columns in the query are to be added in the next release
        # We want the vat creation query to overwrite the destination table because if new data has been put into the pre-vat tables
        # and this workflow has been run an additional time, we dont want duplicates being appended from the original run

        bq --apilog=false query --nouse_legacy_sql --destination_table=~{dataset_name}.~{vat_table} --replace --project_id=~{project_id} \
        'SELECT
            v.vid,
            v.transcript,
            v.contig,
            v.position,
            v.ref_allele,
            v.alt_allele,
            v.gvs_all_ac,
            v.gvs_all_an,
            v.gvs_all_af,
            v.gvs_all_sc,
            v.gvs_max_af,
            v.gvs_max_ac,
            v.gvs_max_an,
            v.gvs_max_sc,
            v.gvs_max_subpop,
            v.gvs_afr_ac,
            v.gvs_afr_an,
            v.gvs_afr_af,
            v.gvs_afr_sc,
            v.gvs_amr_ac,
            v.gvs_amr_an,
            v.gvs_amr_af,
            v.gvs_amr_sc,
            v.gvs_eas_ac,
            v.gvs_eas_an,
            v.gvs_eas_af,
            v.gvs_eas_sc,
            v.gvs_eur_ac,
            v.gvs_eur_an,
            v.gvs_eur_af,
            v.gvs_eur_sc,
            v.gvs_mid_ac,
            v.gvs_mid_an,
            v.gvs_mid_af,
            v.gvs_mid_sc,
            v.gvs_oth_ac,
            v.gvs_oth_an,
            v.gvs_oth_af,
            v.gvs_oth_sc,
            v.gvs_sas_ac,
            v.gvs_sas_an,
            v.gvs_sas_af,
            v.gvs_sas_sc,
            v.gene_symbol,
            v.transcript_source,
            v.aa_change,
            v.consequence,
            v.dna_change_in_transcript,
            v.variant_type,
            v.exon_number,
            v.intron_number,
            v.genomic_location,
            # v.hgvsc AS splice_distance
            v.dbsnp_rsid,
            v.gene_id,
            # v.entrez_gene_id,
            # g.hgnc_gene_id,
            g.gene_omim_id,
            CASE WHEN ( v.transcript is not null and v.is_canonical_transcript is not True)
            THEN False WHEN ( v.transcript is not null and v.is_canonical_transcript is True) THEN True END AS is_canonical_transcript,
            v.gnomad_all_af,
            v.gnomad_all_ac,
            v.gnomad_all_an,
            v.gnomad_failed_filter,
            v.gnomad_max_af,
            v.gnomad_max_ac,
            v.gnomad_max_an,
            v.gnomad_max_subpop,
            v.gnomad_afr_ac,
            v.gnomad_afr_an,
            v.gnomad_afr_af,
            v.gnomad_amr_ac,
            v.gnomad_amr_an,
            v.gnomad_amr_af,
            v.gnomad_asj_ac,
            v.gnomad_asj_an,
            v.gnomad_asj_af,
            v.gnomad_eas_ac,
            v.gnomad_eas_an,
            v.gnomad_eas_af,
            v.gnomad_fin_ac,
            v.gnomad_fin_an,
            v.gnomad_fin_af,
            v.gnomad_nfr_ac,
            v.gnomad_nfr_an,
            v.gnomad_nfr_af,
            v.gnomad_sas_ac,
            v.gnomad_sas_an,
            v.gnomad_sas_af,
            v.gnomad_oth_ac,
            v.gnomad_oth_an,
            v.gnomad_oth_af,
            v.revel,
            v.splice_ai_acceptor_gain_score,
            v.splice_ai_acceptor_gain_distance,
            v.splice_ai_acceptor_loss_score,
            v.splice_ai_acceptor_loss_distance,
            v.splice_ai_donor_gain_score,
            v.splice_ai_donor_gain_distance,
            v.splice_ai_donor_loss_score,
            v.splice_ai_donor_loss_distance,
            g.omim_phenotypes_id,
            g.omim_phenotypes_name,
            v.clinvar_classification,
            v.clinvar_last_updated,
            v.clinvar_phenotype,
        FROM `~{dataset_name}.~{variant_transcript_table}` as v
            left join
        (SELECT gene_symbol, ANY_VALUE(gene_omim_id) AS gene_omim_id, ANY_VALUE(omim_phenotypes_id) AS omim_phenotypes_id, ANY_VALUE(omim_phenotypes_name) AS omim_phenotypes_name FROM `~{dataset_name}.~{genes_table}` group by gene_symbol) as g
        on v.gene_symbol = g.gene_symbol'
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "openbridge/ob_google-bigquery:latest"
        memory: "3 GB"
        preemptible: 3
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        String vat_table_name = vat_table
        Boolean done = true
    }
}
