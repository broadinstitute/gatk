version 1.0

import "GvsExtractCallsetPgen.wdl" as Extract
import "MergePgenHierarchical.wdl" as Merge

workflow GvsExtractCallsetMerged {

    input {
        Boolean go = true
        String dataset_name
        String project_id
        String call_set_identifier

        # Chromosome code corresponding to the set of contig names used by plink that matches those used by the
        # reference
        String pgen_chromosome_code = "chrM"
        # Max number of alt alleles a site can have. If a site exceeds this number, it will not be written
        Int max_alt_alleles = 254
        # If true, does not throw an exception for samples@sites with unsupported ploidy (codes it as missing instead)
        Boolean lenient_ploidy_validation = false

        String cohort_project_id = project_id
        String cohort_dataset_name = dataset_name
        Boolean do_not_filter_override = false
        Boolean control_samples = false
        String extract_table_prefix
        String filter_set_name
        String query_project = project_id
        # This is optional now since the workflow will choose an appropriate value below if this is unspecified.
        Int? scatter_count
        Boolean zero_pad_output_pgen_filenames = true

        # set to "NONE" if all the reference data was loaded into GVS in GvsImportGenomes
        String drop_state = "NONE"

        File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
        Boolean use_interval_weights = true
        File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"

        File? gatk_override
        String? extract_docker_override
        String gatk_docker

        String output_file_base_name = filter_set_name

        Int? extract_maxretries_override
        Int? extract_preemptible_override
        String? output_gcs_dir
        String? split_intervals_extra_args
        Int? split_intervals_disk_size_override
        Int? split_intervals_mem_override
        Float x_bed_weight_scaling = 4
        Float y_bed_weight_scaling = 4
        Boolean write_cost_to_db = true

        # Merge
        String plink_docker
    }

    call Extract.GvsExtractCallset {
        input:
            go = go,
            dataset_name = dataset_name,
            project_id = project_id,
            call_set_identifier = call_set_identifier,
            pgen_chromosome_code = pgen_chromosome_code,
            max_alt_alleles = max_alt_alleles,
            lenient_ploidy_validation = lenient_ploidy_validation,
            cohort_project_id = cohort_project_id,
            cohort_dataset_name = cohort_dataset_name,
            do_not_filter_override = do_not_filter_override,
            control_samples = control_samples,
            extract_table_prefix = extract_table_prefix,
            filter_set_name = filter_set_name,
            query_project = query_project,
            scatter_count = scatter_count,
            zero_pad_output_pgen_filenames = zero_pad_output_pgen_filenames,
            drop_state = drop_state,
            interval_list = interval_list,
            use_interval_weights = use_interval_weights,
            interval_weights_bed = interval_weights_bed,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            output_file_base_name = output_file_base_name,
            extract_maxretries_override = extract_maxretries_override,
            extract_preemptible_override = extract_preemptible_override,
            output_gcs_dir = output_gcs_dir,
            split_intervals_extra_args = split_intervals_extra_args,
            split_intervals_disk_size_override = split_intervals_disk_size_override,
            split_intervals_mem_override = split_intervals_mem_override,
            x_bed_weight_scaling = x_bed_weight_scaling,
            y_bed_weight_scaling = y_bed_weight_scaling,
            write_cost_to_db = write_cost_to_db,
            extract_docker_override = extract_docker_override
    }

    call SplitFilesByChromosome {
        input:
            interval_lists_tar = GvsExtractCallset.output_pgen_interval_files,
            interval_list_filenames = GvsExtractCallset.output_pgen_interval_filenames,
            pgen_files = GvsExtractCallset.output_pgens,
            pvar_files = GvsExtractCallset.output_pvars,
            psam_files = GvsExtractCallset.output_psams,
            split_intervals_disk_size_override = split_intervals_disk_size_override
    }

    scatter(i in range(length(SplitFilesByChromosome.pgen_lists))) {
        Int split_count = ceil(length(GvsExtractCallset.output_pgens)/(length(SplitFilesByChromosome.pgen_lists)*10))
        String contig = sub(basename(SplitFilesByChromosome.pgen_lists[i]), ".pgen_list", "")
        call Merge.MergePgenWorkflow {
            input:
                pgen_file_list = SplitFilesByChromosome.pgen_lists[i],
                pvar_file_list = SplitFilesByChromosome.pvar_lists[i],
                psam_file_list = SplitFilesByChromosome.psam_lists[i],
                plink_docker = plink_docker,
                output_file_base_name = "~{output_file_base_name}.${contig}",
                merge_disk_size = 1024,
                split_count = split_count,
                zero_padded_prefix = zero_pad_output_pgen_filenames
        }
    }

    output {
        Array[File] merged_pgens = MergePgenWorkflow.pgen_file
        Array[File] merged_pvars = MergePgenWorkflow.pvar_file
        Array[File] merged_psams = MergePgenWorkflow.psam_file
        Array[File] output_pgens = GvsExtractCallset.output_pgens
        Array[File] output_pvars = GvsExtractCallset.output_pvars
        Array[File] output_psams = GvsExtractCallset.output_psams
        File output_pgen_interval_files = GvsExtractCallset.output_pgen_interval_files
        Float total_pgens_size_mb = GvsExtractCallset.total_pgens_size_mb
        File manifest = GvsExtractCallset.manifest
        File? sample_name_list = GvsExtractCallset.sample_name_list
    }

}

task SplitFilesByChromosome {

    input {
        File interval_lists_tar
        Array[String] interval_list_filenames
        Array[String] pgen_files
        Array[String] pvar_files
        Array[String] psam_files

        Int? split_intervals_disk_size_override
    }

    Int disk_size = select_first([split_intervals_disk_size_override, 10])

    command <<<
        set -euxo pipefail

        PGEN_ARRAY=(~{sep=" " pgen_files})
        PSAM_ARRAY=(~{sep=" " psam_files})
        PVAR_ARRAY=(~{sep=" " pvar_files})

        INTERVAL_LIST_ARRAY=(~{sep=" " interval_list_filenames})

        # Extract the interval lists from the tar file
        tar -xvf ~{interval_lists_tar}

        # Loop through the interval lists, extracting the chromosome code for each, and writing the
        # corresponding pgen, pvar, and psam files to a separate list file for each chromosome
        for index in "${!INTERVAL_LIST_ARRAY[@]}"; do
            # Extract the chromosome code from the interval list by getting from the first interval in the file
            chrom=$(awk '!/^@/ {print $1; exit}' ${INTERVAL_LIST_ARRAY[$index]})
            # Append the corresponding pgen, pvar, and psam files to a list file for this chromosome
            echo "${PGEN_ARRAY[$index]}" >> ${chrom}.pgen_list
            echo "${PVAR_ARRAY[$index]}" >> ${chrom}.pvar_list
            echo "${PSAM_ARRAY[$index]}" >> ${chrom}.psam_list
        done
    >>>

    output {
        Array[File] pgen_lists = glob("*.pgen_list")
        Array[File] pvar_lists = glob("*.pvar_list")
        Array[File] psam_lists = glob("*.psam_list")
    }
    
    runtime {
        docker: "ubuntu:20.04"
        memory: "1 GB"
        disks: "local-disk ${disk_size} HDD"
        cpu: "1"
        preemptible: 1
    }
}