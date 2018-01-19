# interval list .tsv file column names
contig_column_name = 'CONTIG'
start_column_name = 'START'
end_column_name = 'END'
count_column_name = 'COUNT'

# prefix for saving posteriors for multiple samples
sample_folder_prefix = "SAMPLE_"

# log copy number posterior matrix column name prefix for each integer copy number state
copy_number_column_prefix = "COPY_NUMBER_"

# ploidy prior table header column names
ploidy_prior_contig_name_column = 'CONTIG_NAME'
ploidy_prior_prefix = 'PLOIDY_PRIOR_'

# column names for ploidy and depth .tsv outputs
sample_name_column_name = 'SAMPLE_NAME'
global_read_depth_column_name = 'GLOBAL_READ_DEPTH'
average_ploidy_column_name = 'AVERAGE_PLOIDY'
ploidy_column_name = 'PLOIDY'
ploidy_gq_column_name = 'PLOIDY_GQ'

# column names for copy-number segments file
copy_number_call_column_name = 'COPY_NUMBER_CALL'
num_spanning_intervals_column_name = 'NUM_SPANNING_INTERVALS'
some_quality_column_name = 'SOME_QUALITY'
exact_quality_column_name = 'EXACT_QUALITY'
start_quality_column_name = 'START_QUALITY'
end_quality_column_name =  'END_QUALITY'

# regular expression for matching sample name from header comment line
sample_name_header_regexp = "^@RG.*SM:(.*)[\t]*.*$"

# prefix for adding sample name as a header comment line
sample_name_sam_header_prefix = "RG\tID:GATKCopyNumber\tSM:"

# default file names for loading and saving models, posteriors, and configurations
default_sample_read_depth_tsv_filename = 'global_read_depth.tsv'
default_sample_name_txt_filename = 'sample_name.txt'
default_sample_contig_ploidy_tsv_filename = 'contig_ploidy.tsv'
default_copy_number_log_posterior_tsv_filename = "log_q_c_tc.tsv"
default_copy_number_log_emission_tsv_filename = "log_c_emission_tc.tsv"
default_class_log_posterior_tsv_filename = "log_q_tau_tk.tsv"
default_baseline_copy_number_tsv_filename = "baseline_copy_number_t.tsv"
default_copy_number_segments_tsv_filename = "copy_number_segments.tsv"

default_denoising_config_json_filename = "denoising_config.json"
default_calling_config_json_filename = "calling_config.json"
default_ploidy_config_json_filename = "ploidy_config.json"
default_gcnvkernel_version_json_filename = "gcnvkernel_version.json"

default_interval_list_filename = "interval_list.tsv"
default_contig_ploidy_prior_tsv_filename = 'contig_ploidy_prior.tsv'

default_adamax_m_filename = "adamax_m.npy"
default_adamax_u_filename = "adamax_u.npy"
default_adamax_res_filename = "adamax_res.npy"
