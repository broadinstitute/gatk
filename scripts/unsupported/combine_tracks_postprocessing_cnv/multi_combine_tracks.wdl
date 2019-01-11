# An unsupported workflow for post-processing CNV calls and making the files IGV- and absolute-compatible.
#  Additionally, this will also create one seg file (compatible with IGV) for all samples specified.
import "combine_tracks.wdl" as CombineTracks
import "aggregate_combined_tracks.wdl" as AggregateCombineTracks

workflow MultiCombineTracks {
    	Array[File] tumor_called_segs
    	Array[File] matched_normal_called_segs
    	Array[File] tumor_modeled_segs
        Array[File] af_params
        File ref_fasta
        File ref_fasta_dict
        File ref_fasta_fai
        File centromere_tracks_seg
        File gistic_blacklist_tracks_seg
        File? gatk4_jar_override
        Array[String] columns_of_interest
        Int? germline_tagging_padding
        String group_id
        String gatk_docker
        Int? min_hets_acs_results

        scatter (i in range(length(tumor_called_segs))) {

            call CombineTracks.CombineTracksWorkflow as CombineTracksWorkflow {
                input:
                	tumor_called_seg = tumor_called_segs[i],
                	tumor_modeled_seg = tumor_modeled_segs[i],
                	af_param = af_params[i],
                	matched_normal_called_seg = matched_normal_called_segs[i],
                    ref_fasta = ref_fasta,
                    ref_fasta_dict = ref_fasta_dict,
                    ref_fasta_fai = ref_fasta_fai,
                    centromere_tracks_seg = centromere_tracks_seg,
                    gistic_blacklist_tracks_seg = gistic_blacklist_tracks_seg,
                    gatk4_jar_override = gatk4_jar_override,
                    columns_of_interest = columns_of_interest,
                    germline_tagging_padding = germline_tagging_padding,
                    group_id = group_id,
                    gatk_docker = gatk_docker,
                    min_hets_acs_results = min_hets_acs_results
            }
        }

        call AggregateCombineTracks.AggregateCombinedTracksWorkflow as Aggregate {
            input:
                    group_id = group_id,
                    tumor_with_germline_filtered_segs = CombineTracksWorkflow.cnv_postprocessing_tumor_with_tracks_filtered_merged_seg,
                    normals_igv_compat = CombineTracksWorkflow.cnv_postprocessing_normal_igv_compat,
                    tumors_igv_compat = CombineTracksWorkflow.cnv_postprocessing_tumor_igv_compat,
                    tumors_gistic2_compat = CombineTracksWorkflow.cnv_postprocessing_tumor_with_tracks_filtered_merged_seg_gistic2
        }
        output {
            File tumor_with_germline_filtered_segs = Aggregate.cnv_postprocessing_aggregated_tumors_post
            File normals_igv_compat = Aggregate.cnv_postprocessing_aggregated_normals
            File tumors_igv_compat = Aggregate.cnv_postprocessing_aggregated_tumors_pre
            File tumor_with_germline_filtered_segs_gistic2 = Aggregate.cnv_postprocessing_aggregated_tumors_post_gistic2
            Array[File] tumor_acs_compat = CombineTracksWorkflow.cnv_postprocessing_tumor_acs_seg
            Array[File] tumor_acs_skew = CombineTracksWorkflow.cnv_postprocessing_tumor_acs_skew
        }
}