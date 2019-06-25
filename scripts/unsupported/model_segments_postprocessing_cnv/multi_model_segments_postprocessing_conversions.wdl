# TODO: For docs, assumes that the blacklist was applied in the running of ModelSegments (`-XL`)
import "model_segments_postprocessing_conversions.wdl" as ModelSegmentsPostProcessingConversions
#import "aggregate_model_segments_postprocessing.wdl" as AggregateModelSegmentsPostProcessing

workflow MultiModelSegmentsPostProcessingConversions {
    	Array[File] tumor_called_segs
    	Array[File] tumor_modeled_segs
        Array[File] af_params
        File ref_fasta
        File ref_fasta_dict
        File ref_fasta_fai
        String group_id
        Float? maf90_threshold
        Float maf90_threshold_final = select_first([maf90_threshold, 0.485])
        Int? min_hets_acs_results
        File? gatk4_jar_override
        String gatk_docker
        Array[String]? columns_of_interest

        scatter (i in range(length(tumor_called_segs))) {

            call ModelSegmentsPostProcessingConversions.ModelSegmentsPostProcessingConversionWorkflow as ModelSegmentsPostProcessingConversionWorkflow {
                input:
                	tumor_called_seg = tumor_called_segs[i],
                	tumor_modeled_seg = tumor_modeled_segs[i],
                    ref_fasta = ref_fasta,
                    ref_fasta_dict = ref_fasta_dict,
                    ref_fasta_fai = ref_fasta_fai,
                    af_param = af_params[i],
                    min_hets_acs_results = min_hets_acs_results,
                    maf90_threshold = maf90_threshold_final,
                    gatk_docker = gatk_docker,
                    columns_of_interest = columns_of_interest
            }
        }

#        call AggregateModelSegmentsPostProcessing.AggregateModelSegmentsPostProcessingWorkflow as Aggregate {
#            input:
#                    group_id = group_id,
#                    tumor_with_germline_filtered_segs = ModelSegmentsPostProcessingConversionWorkflow.cnv_postprocessing_tumor_with_tracks_filtered_merged_seg,
#                    normals_igv_compat = ModelSegmentsPostProcessingConversionWorkflow.cnv_postprocessing_normal_igv_compat,
#                    tumors_igv_compat = ModelSegmentsPostProcessingConversionWorkflow.cnv_postprocessing_tumor_igv_compat,
#                    tumors_gistic2_compat = ModelSegmentsPostProcessingConversionWorkflow.cnv_postprocessing_tumor_with_tracks_filtered_merged_seg_gistic2
#        }
        output {
#            File tumor_with_germline_filtered_segs = Aggregate.cnv_postprocessing_aggregated_tumors_post
#            File normals_igv_compat = Aggregate.cnv_postprocessing_aggregated_normals
#            File tumors_igv_compat = Aggregate.cnv_postprocessing_aggregated_tumors_pre
#            File tumor_with_germline_filtered_segs_gistic2 = Aggregate.cnv_postprocessing_aggregated_tumors_post_gistic2
            Array[File] tumor_acs_compat = ModelSegmentsPostProcessingConversionWorkflow.cnv_postprocessing_tumor_acs_seg
            Array[File] tumor_acs_skew = ModelSegmentsPostProcessingConversionWorkflow.cnv_postprocessing_tumor_acs_skew
        }
}