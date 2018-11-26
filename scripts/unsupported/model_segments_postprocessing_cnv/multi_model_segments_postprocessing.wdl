# An unsupported workflow for post-processing CNV calls and making the files IGV- and absolute-compatible.
#  Additionally, this will also create one seg file (compatible with IGV) for all samples specified.
import "model_segments_postprocessing.wdl" as ModelSegmentsPostProcessing
import "aggregate_model_segments_postprocessing.wdl" as AggregateModelSegmentsPostProcessing

workflow MultiModelSegmentsPostProcessing {
    	Array[File] tumor_called_segs
    	Array[File] matched_normal_called_segs
    	Array[File] tumor_modeled_segs
    	Array[File] matched_normal_modeled_segs
        Array[File] af_params
        File ref_fasta
        File ref_fasta_dict
        File ref_fasta_fai
        File centromere_tracks_seg
        File gistic_blacklist_tracks_seg
        File? gatk4_jar_override
        Int? germline_tagging_padding
        String group_id
        String gatk_docker
        Int? min_hets_acs_results
        Boolean? is_ignore_cnloh_in_matched_normal
        Int? cnloh_check_max_size
        scatter (i in range(length(tumor_called_segs))) {

            call ModelSegmentsPostProcessing.ModelSegmentsPostProcessingWorkflow as ModelSegmentsPostProcessingWorkflow {
                input:
                	tumor_called_seg = tumor_called_segs[i],
                	tumor_modeled_seg = tumor_modeled_segs[i],
                	af_param = af_params[i],
                	matched_normal_called_seg = matched_normal_called_segs[i],
                	matched_normal_modeled_seg = matched_normal_modeled_segs[i],
                    ref_fasta = ref_fasta,
                    ref_fasta_dict = ref_fasta_dict,
                    ref_fasta_fai = ref_fasta_fai,
                    centromere_tracks_seg = centromere_tracks_seg,
                    gistic_blacklist_tracks_seg = gistic_blacklist_tracks_seg,
                    gatk4_jar_override = gatk4_jar_override,
                    germline_tagging_padding = germline_tagging_padding,
                    group_id = group_id,
                    gatk_docker = gatk_docker,
                    min_hets_acs_results = min_hets_acs_results,
                    is_ignore_cnloh_in_matched_normal = is_ignore_cnloh_in_matched_normal,
                    cnloh_check_max_size = cnloh_check_max_size
            }
        }

        call AggregateModelSegmentsPostProcessing.AggregateModelSegmentsPostProcessingWorkflow as Aggregate {
            input:
                    group_id = group_id,
                    tumor_with_germline_filtered_segs = ModelSegmentsPostProcessingWorkflow.cnv_postprocessing_tumor_with_tracks_filtered_merged_seg,
                    normals_igv_compat = ModelSegmentsPostProcessingWorkflow.cnv_postprocessing_normal_igv_compat,
                    tumors_igv_compat = ModelSegmentsPostProcessingWorkflow.cnv_postprocessing_tumor_igv_compat,
                    tumors_gistic2_compat = ModelSegmentsPostProcessingWorkflow.cnv_postprocessing_tumor_with_tracks_filtered_merged_seg_gistic2
        }
        output {
            File tumor_with_germline_filtered_segs = Aggregate.cnv_postprocessing_aggregated_tumors_post
            File normals_igv_compat = Aggregate.cnv_postprocessing_aggregated_normals
            File tumors_igv_compat = Aggregate.cnv_postprocessing_aggregated_tumors_pre
            File tumor_with_germline_filtered_segs_gistic2 = Aggregate.cnv_postprocessing_aggregated_tumors_post_gistic2
            Array[File] tumor_acs_compat = ModelSegmentsPostProcessingWorkflow.cnv_postprocessing_tumor_acs_seg
            Array[File] tumor_acs_skew = ModelSegmentsPostProcessingWorkflow.cnv_postprocessing_tumor_acs_skew
        }
}