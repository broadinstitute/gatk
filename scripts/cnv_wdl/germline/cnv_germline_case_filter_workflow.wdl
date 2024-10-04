version 1.0

import "cnv_germline_case_workflow.wdl" as GermlineCNVCaseWorkflow

workflow SingleSampleGCNVAndFilterVCFs {
    input {
        Array[String]? filter_expressions #= "QUAL < 50"  # Example filter criteria
        Array[String]? filter_names

        ##################################
        ##for case mode
        #### required basic arguments ####
        ##################################
        File intervals
        File? blacklist_intervals
        File filtered_intervals
        String normal_bam
        String normal_bai
        File contig_ploidy_model_tar
        Array[File]+ gcnv_model_tars
        Array[File]+ pon_genotyped_segments_vcfs
        Int num_intervals_per_scatter
        File ref_fasta_dict
        File ref_fasta_fai
        File ref_fasta
        String gatk_docker

        ##################################
        #### optional basic arguments ####
        ##################################
        File? gatk4_jar_override
        Int? preemptible_attempts

        # Required if BAM/CRAM is in a requester pays bucket
        String? gcs_project_for_requester_pays

        ####################################################
        #### optional arguments for PreprocessIntervals ####
        ####################################################
        Int? padding
        Int? bin_length

        ##############################################
        #### optional arguments for CollectCounts ####
        ##############################################
        Array[String]? disabled_read_filters_for_collect_counts
        String? collect_counts_format
        Boolean? collect_counts_enable_indexing
        Int? mem_gb_for_collect_counts

        ######################################################################
        #### optional arguments for DetermineGermlineContigPloidyCaseMode ####
        ######################################################################
        Float? ploidy_mapping_error_rate
        Float? ploidy_sample_psi_scale
        Int? mem_gb_for_determine_germline_contig_ploidy
        Int? cpu_for_determine_germline_contig_ploidy
        Int? disk_for_determine_germline_contig_ploidy

        ##########################################################
        #### optional arguments for GermlineCNVCallerCaseMode ####
        ##########################################################
        Float? gcnv_p_alt
        Float? gcnv_cnv_coherence_length
        Int? gcnv_max_copy_number
        Int? mem_gb_for_germline_cnv_caller
        Int? cpu_for_germline_cnv_caller
        Int? disk_for_germline_cnv_caller

        # optional arguments for germline CNV denoising model
        Float? gcnv_mapping_error_rate
        Float? gcnv_sample_psi_scale
        Float? gcnv_depth_correction_tau
        String? gcnv_copy_number_posterior_expectation_mode
        Int? gcnv_active_class_padding_hybrid_mode

        # optional arguments for Hybrid ADVI
        Float? gcnv_learning_rate
        Float? gcnv_adamax_beta_1
        Float? gcnv_adamax_beta_2
        Int? gcnv_log_emission_samples_per_round
        Float? gcnv_log_emission_sampling_median_rel_error
        Int? gcnv_log_emission_sampling_rounds
        Int? gcnv_max_advi_iter_first_epoch
        Int? gcnv_max_advi_iter_subsequent_epochs
        Int? gcnv_min_training_epochs
        Int? gcnv_max_training_epochs
        Float? gcnv_initial_temperature
        Int? gcnv_num_thermal_advi_iters
        Int? gcnv_convergence_snr_averaging_window
        Float? gcnv_convergence_snr_trigger_threshold
        Int? gcnv_convergence_snr_countdown_window
        Int? gcnv_max_calling_iters
        Float? gcnv_caller_update_convergence_threshold
        Float? gcnv_caller_internal_admixing_rate
        Float? gcnv_caller_external_admixing_rate
        Boolean? gcnv_disable_annealing

        ###################################################
        #### arguments for PostprocessGermlineCNVCalls ####
        ###################################################
        Int ref_copy_number_autosomal_contigs
        Array[String]? allosomal_contigs

        ##########################
        #### arguments for QC ####
        ##########################
        Int maximum_number_events_per_sample
        Int maximum_number_pass_events_per_sample

    }


    call GermlineCNVCaseWorkflow.CNVGermlineCaseWorkflow as CNVGermlineCaseWorkflow {
        input:
            intervals = intervals,
            blacklist_intervals=blacklist_intervals,
            filtered_intervals=filtered_intervals,
            normal_bams=[normal_bam],
            normal_bais=[normal_bai],
            contig_ploidy_model_tar=contig_ploidy_model_tar,
            gcnv_model_tars=gcnv_model_tars,
            num_intervals_per_scatter=num_intervals_per_scatter,
            ref_fasta_dict=ref_fasta_dict,
            ref_fasta_fai=ref_fasta_fai,
            ref_fasta=ref_fasta,
            gatk_docker=gatk_docker,

            ##################################
            #### optional basic arguments ####
            ##################################
            gatk4_jar_override=gatk4_jar_override,
            preemptible_attempts=preemptible_attempts,

            # Required if BAM/CRAM is in a requester pays bucket
            gcs_project_for_requester_pays=gcs_project_for_requester_pays,

            ####################################################
            #### optional arguments for PreprocessIntervals ####
            ####################################################
            padding=padding,
            bin_length=bin_length,

            ##############################################
            #### optional arguments for CollectCounts ####
            ##############################################
            disabled_read_filters_for_collect_counts=disabled_read_filters_for_collect_counts,
            collect_counts_format=collect_counts_format,
            collect_counts_enable_indexing=collect_counts_enable_indexing,
            mem_gb_for_collect_counts=mem_gb_for_collect_counts,

            ######################################################################
            #### optional arguments for DetermineGermlineContigPloidyCaseMode ####
            ######################################################################
            ploidy_mapping_error_rate=ploidy_mapping_error_rate,
            ploidy_sample_psi_scale=ploidy_sample_psi_scale,
            mem_gb_for_determine_germline_contig_ploidy=mem_gb_for_determine_germline_contig_ploidy,
            cpu_for_determine_germline_contig_ploidy=cpu_for_determine_germline_contig_ploidy,
            disk_for_determine_germline_contig_ploidy=disk_for_determine_germline_contig_ploidy,

            ##########################################################
            #### optional arguments for GermlineCNVCallerCaseMode ####
            ##########################################################
            gcnv_p_alt=gcnv_p_alt,
            gcnv_cnv_coherence_length=gcnv_cnv_coherence_length,
            gcnv_max_copy_number=gcnv_max_copy_number,
            mem_gb_for_germline_cnv_caller=mem_gb_for_germline_cnv_caller,
            cpu_for_germline_cnv_caller=cpu_for_germline_cnv_caller,
            disk_for_germline_cnv_caller=disk_for_germline_cnv_caller,

            # optional arguments for germline CNV denoising model
            gcnv_mapping_error_rate=gcnv_mapping_error_rate,
            gcnv_sample_psi_scale=gcnv_sample_psi_scale,
            gcnv_depth_correction_tau=gcnv_depth_correction_tau,
            gcnv_copy_number_posterior_expectation_mode=gcnv_copy_number_posterior_expectation_mode,
            gcnv_active_class_padding_hybrid_mode=gcnv_active_class_padding_hybrid_mode,

            # optional arguments for Hybrid ADVI
            gcnv_learning_rate=gcnv_learning_rate,
            gcnv_adamax_beta_1=gcnv_adamax_beta_1,
            gcnv_adamax_beta_2=gcnv_adamax_beta_2,
            gcnv_log_emission_samples_per_round=gcnv_log_emission_samples_per_round,
            gcnv_log_emission_sampling_median_rel_error=gcnv_log_emission_sampling_median_rel_error,
            gcnv_log_emission_sampling_rounds=gcnv_log_emission_sampling_rounds,
            gcnv_max_advi_iter_first_epoch=gcnv_max_advi_iter_first_epoch,
            gcnv_max_advi_iter_subsequent_epochs=gcnv_max_advi_iter_subsequent_epochs,
            gcnv_min_training_epochs=gcnv_min_training_epochs,
            gcnv_max_training_epochs=gcnv_max_training_epochs,
            gcnv_initial_temperature=gcnv_initial_temperature,
            gcnv_num_thermal_advi_iters=gcnv_num_thermal_advi_iters,
            gcnv_convergence_snr_averaging_window=gcnv_convergence_snr_averaging_window,
            gcnv_convergence_snr_trigger_threshold=gcnv_convergence_snr_trigger_threshold,
            gcnv_convergence_snr_countdown_window=gcnv_convergence_snr_countdown_window,
            gcnv_max_calling_iters=gcnv_max_calling_iters,
            gcnv_caller_update_convergence_threshold=gcnv_caller_update_convergence_threshold,
            gcnv_caller_internal_admixing_rate=gcnv_caller_internal_admixing_rate,
            gcnv_caller_external_admixing_rate=gcnv_caller_external_admixing_rate,
            gcnv_disable_annealing=gcnv_disable_annealing,

            ###################################################
            #### arguments for PostprocessGermlineCNVCalls ####
            ###################################################
            ref_copy_number_autosomal_contigs=ref_copy_number_autosomal_contigs,
            allosomal_contigs=allosomal_contigs,

            ##########################
            #### arguments for QC ####
            ##########################
            maximum_number_events_per_sample=maximum_number_events_per_sample,
            maximum_number_pass_events_per_sample=maximum_number_pass_events_per_sample
    }

    File vcf = CNVGermlineCaseWorkflow.genotyped_segments_vcfs[0]
    File vcf_index = CNVGermlineCaseWorkflow.genotyped_intervals_vcf_indexes[0]
    call ExtractSamplename {
        input:
            vcf_filename = basename(vcf),
            gatk_docker = gatk_docker
    }

    call ExtractPoNFreq {
        input:
            vcf = vcf,
            panel_vcfs = pon_genotyped_segments_vcfs,
            intervals = intervals
    }

    call AnnotateWithPoNFreq {
        input:
            vcf = vcf,
            vcf_idx = vcf_index,
            annotations = ExtractPoNFreq.annotations
    }

    call FilterVCF{
        input:
            samplename = ExtractSamplename.samplename,
            vcf_file = AnnotateWithPoNFreq.output_vcf,
            filter_expressions = filter_expressions,
            filter_names = filter_names,
            gatk_docker = gatk_docker
    }



    output {
        File filtered_vcf = FilterVCF.filtered_vcf
        File filtered_vcf_index = FilterVCF.filtered_vcf
        File filtered_vcf_md5sum = FilterVCF.filtered_vcf_md5sum

        File preprocessed_intervals = CNVGermlineCaseWorkflow.preprocessed_intervals
        File read_counts_entity_id = CNVGermlineCaseWorkflow.read_counts_entity_id[0]
        File read_counts = CNVGermlineCaseWorkflow.read_counts[0]
        File sample_contig_ploidy_calls_tars = CNVGermlineCaseWorkflow.sample_contig_ploidy_calls_tars[0]
        Array[File] gcnv_calls_tars = CNVGermlineCaseWorkflow.gcnv_calls_tars[0]
        File gcnv_tracking_tars = CNVGermlineCaseWorkflow.gcnv_tracking_tars[0]
        File genotyped_intervals_vcfs = CNVGermlineCaseWorkflow.genotyped_intervals_vcfs[0]
        File genotyped_intervals_vcf_indexes = CNVGermlineCaseWorkflow.genotyped_intervals_vcf_indexes[0]
        File genotyped_segments_vcfs = CNVGermlineCaseWorkflow.genotyped_segments_vcfs[0]
        File genotyped_segments_vcf_indexes = CNVGermlineCaseWorkflow.genotyped_segments_vcf_indexes[0]
        File qc_status_files = CNVGermlineCaseWorkflow.qc_status_files[0]
        String qc_status_strings = CNVGermlineCaseWorkflow.qc_status_strings[0]
        File denoised_copy_ratios = CNVGermlineCaseWorkflow.denoised_copy_ratios[0]
    }
}


task ExtractSamplename {
    input {
        String vcf_filename
        String gatk_docker
    }
    command <<<
        echo ~{vcf_filename} | sed -E 's/genotyped-segments-(.*)\.cram\.vcf\.gz/\1/' > output.txt
        >>>
    output {
        String samplename = read_string("output.txt")
    }
    runtime {
        docker: gatk_docker
        memory: "2G"
        cpu: 2
    }
}

task ExtractPoNFreq {
    input {
        File vcf
        Array[File] panel_vcfs
        File intervals
        Float overlap_thresh = 0.5
        Int mem_gb = 4
        Int disk_size_gb = 100
    }

    String basename_out = basename(vcf, ".vcf.gz")

    command <<<
        python << "EOF"
        import pandas as pd
        import numpy as np

        intervals = pd.read_csv("~{intervals}", sep="\t", comment="@", names = ["contig","start","end","dummy1","dummy2"])

        def add_exon_idxs(df, exons):
            contigs = set(df.contig)
            for contig in contigs:
                df.loc[df.contig==contig,"start_exon_idx"]=np.searchsorted(exons.loc[exons.contig==contig].start,
                                                                           df.loc[df.contig==contig].start,"left")
                df.loc[df.contig==contig,"end_exon_idx"]=np.searchsorted(exons.loc[exons.contig==contig].start,
                                                                           df.loc[df.contig==contig].end,"right")

        def standardize_gt_vcf(df):
            df['end'] = df['INFO'].str.replace('END=','').astype(int) #split('=').apply(lambda x: x[1])
            df["svtype"] = df.ALT.str.replace("<","").str.replace(">","")
            for num,f in enumerate(df.loc[0,'FORMAT'].split(':')):
                df[f] = df['SAMPLE'].str.split(':').apply(lambda x: x[num])
            return df
        
        def read_vcf_to_df(vcf):
            df = pd.read_csv(vcf, sep='\t',comment='#',compression='gzip',
                                 names=['contig','start','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE'])
            df = standardize_gt_vcf(df)
            df['sample_path']=vcf
            for c in ['NP','CN','QA','QS']:
                df[c]=df[c].astype(int)
            return df
            
        df = read_vcf_to_df("~{vcf}")
        add_exon_idxs(df, intervals)
        df = df.sort_values(["contig","start_exon_idx","end_exon_idx"])

        vcfs = ["~{sep='","' panel_vcfs}"]
        df_panel = pd.concat([read_vcf_to_df(vcf) for vcf in vcfs])
        add_exon_idxs(df_panel, intervals)
        df_panel = df_panel.sort_values(["contig","start_exon_idx","end_exon_idx"])

        def annotate_with_panel_count(panel_df, df, svtype, contig):
            panel_df_mask = ((panel_df.svtype==svtype) &
                             (panel_df.contig == contig))
            df_mask = ((df.svtype==svtype) &
                             (df.contig == contig)
                          )
            panel_df_slice = panel_df.loc[panel_df_mask]
            df_slice = df.loc[df_mask]

            df.loc[df_mask,"start_panel_idx"]=np.searchsorted(panel_df_slice.end_exon_idx,
                                                                 df_slice.start_exon_idx, "left")
            df.loc[df_mask,"end_panel_idx"]=np.searchsorted(panel_df_slice.start_exon_idx,
                                                                 df_slice.end_exon_idx, "right")

            df.loc[df_mask & (df.start_panel_idx ==
                                         df.end_panel_idx), "PANEL_COUNT"] = 0
            panel_counts=[]
            for i, r in df.loc[df_mask & (df.start_panel_idx !=
                                         df.end_panel_idx)].iterrows():
                overlapping_panel_events = panel_df_slice.iloc[int(r.start_panel_idx):int(r.end_panel_idx)]
                overlapping_exon_lengths = (np.minimum(r.end_exon_idx,overlapping_panel_events.end_exon_idx) -
                                           np.maximum(r.start_exon_idx,overlapping_panel_events.start_exon_idx))
                overlapping_panel_fracs = overlapping_exon_lengths/(r.end_exon_idx - r.start_exon_idx)
                overlapping_panel_events = overlapping_panel_events.loc[overlapping_panel_fracs>=0.5]
                panel_count =  len(set(overlapping_panel_events.sample_path))
                panel_counts.append(panel_count)

            df.loc[df_mask & (df.start_panel_idx !=
                                         df.end_panel_idx), "PANEL_COUNT"] = panel_counts

        n_panel_samples = len(set(df_panel.sample_path))
        contigs = set(df.contig)
        svtypes = set(df.svtype)
        for contig in contigs:
            for svtype in svtypes:
                annotate_with_panel_count(df_panel, df, svtype, contig)

        df.loc[:,"PANEL_FREQ"]=(df.PANEL_COUNT/n_panel_samples)

        df_annotations = df[["contig","start","PANEL_FREQ","PANEL_COUNT"]].copy()
        df_annotations = df_annotations.rename({"contig":"CHROM","start":"POS"}, axis=1)
        df_annotations = df_annotations.astype({"PANEL_COUNT":int})
        df_annotations.to_csv("~{basename_out}.annotations.tsv", index = False, sep="\t", float_format="%g")


        EOF

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim"
        preemptible: 3
        cpu: 2
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GB"
    }

    output {
        File annotations = "~{basename_out}.annotations.tsv"
    }
}

task AnnotateWithPoNFreq {
    input {
        File vcf
        File vcf_idx
        File annotations
        Int mem_gb=4
        Int disk_size_gb = 100
    }

    String output_basename = basename(vcf)
    command <<<
        bgzip ~{annotations}
        tabix -s1 -b2 -e2 --skip-lines 1 ~{annotations}.gz

        bcftools annotate -a ~{annotations}.gz -c CHROM,POS,PANEL_FREQ,PANEL_COUNT -o  ~{output_basename}.vcf.gz ~{vcf}

    >>>

    runtime {
            docker: "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
            disks: "local-disk " + disk_size_gb + " HDD"
            memory: mem_gb + " GiB"
    }

    output {
        File output_vcf = "~{output_basename}.vcf.gz"
    }
}


task FilterVCF {
    input {
        String samplename
        File vcf_file
        Array[String] filter_expressions = ["(FMT/CN>2 & QUAL<50) | (FMT/CN==1 & QUAL<100 ) | (FMT/CN==0 & QUAL<400)","INFO/PANEL_COUNT>1"]
        String gatk_docker
        Array[String] filter_names = ["LowQual","PanelCount"]
    }

    command <<<
        cp ~{vcf_file} tmp.vcf.gz
        filters=('~{sep="' '" filter_expressions}')
        fitler_names=(~{sep=" " filter_names}

        for i in ${!filters[@]}
        do
            eval bcftools filter -m + -e \'${filters[$i]}\' --soft-filter ${filter_names[$i]} -Oz -o tmp_out.vcf.gz tmp.vcf.gz
            mv tmp_out.vcf.gz tmp.vcf.gz
        done

        mv tmp.vcf.gz ~{samplename}.filtered.genotyped-segments.vcf.gz

        bcftools index -t ~{samplename}.filtered.genotyped-segments.vcf.gz

        md5sum ~{samplename}.filtered.genotyped-segments.vcf.gz | awk '{ print $1 }' > ~{samplename}.filtered.genotyped-segments.vcf.gz.md5sum
    >>>

    output {
        File filtered_vcf = samplename + ".filtered.genotyped-segments.vcf.gz"
        File filtered_vcf_index = samplename + ".filtered.genotyped-segments.vcf.gz.tbi"
        File filtered_vcf_md5sum = samplename + ".filtered.genotyped-segments.vcf.gz.md5sum"
    }

    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        memory: "8G"
        cpu: 2
    }
}