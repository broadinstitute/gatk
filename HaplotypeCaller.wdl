version 1.0

# Run HaplotypeCaller (GATK Version 4.1.5.0-27-g4c9f5fa-SNAPSHOT)
#
# Call germline SNPs and indels via local re-assembly of haplotypes
#
#  General Workflow (non-tool) Arguments
#    gatk                                              Location of gatk used to run this workflow
#
#  Required Tool Arguments
#    I                                                  BAM/SAM/CRAM file containing reads                          
#    O                                                  File to which variants should be written                    
#    reference                                          Reference sequence file                                     
#
#  Optional Tool Arguments
#    alleles                                            The set of alleles to force-call regardless of evidence     
#    annotate_with_num_discovered_alleles               If provided, we will annotate records with the number of alternate alleles that 
#    annotation                                         One or more specific annotations to add to variant calls    
#    annotation_group                                   One or more groups of annotations to apply to variant calls 
#    annotations_to_exclude                             One or more specific annotations to exclude from variant calls
#    arguments_file                                     read one or more arguments files and add them to the command line
#    assembly_region_out                                Output the assembly region to this IGV formatted file       
#    assembly_region_padding                            Number of additional bases of context to include around each assembly region
#    base_quality_score_threshold                       Base qualities below this threshold will be reduced to the minimum (6)
#    cloud_index_prefetch_buffer                        Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudP
#    cloud_prefetch_buffer                              Size of the cloud-only prefetch buffer (in MB; 0 to disable).
#    contamination_fraction_to_filter                   Fraction of contamination in sequencing data (for all samples) to aggressively r
#    correct_overlapping_quality                        Undocumented option                                         
#    dbsnp                                              dbSNP file                                                  
#    disable_bam_index_caching                          If true, don't cache bam indexes, this will reduce memory requirements but may h
#    disable_sequence_dictionary_validation             If specified, do not check the sequence dictionaries from our inputs for compati
#    founder_id                                         Samples representing the population "founders"              
#    gcs_max_retries                                    If the GCS bucket channel errors out, how many times it will attempt to re-initi
#    gcs_project_for_requester_pays                     Project to bill when accessing "requester pays" buckets. If unset, these buckets
#    graph_output                                       Write debug assembly graph information to this file         
#    help                                               display the help message                                    
#    heterozygosity                                     Heterozygosity value used to compute prior likelihoods for any locus.  See the G
#    heterozygosity_stdev                               Standard deviation of heterozygosity for SNP and indel calling.
#    indel_heterozygosity                               Heterozygosity for indel calling.  See the GATKDocs for heterozygosity for full 
#    interval_merging_rule                              Interval merging rule for abutting intervals                
#    intervals                                          One or more genomic intervals over which to operate         
#    max_assembly_region_size                           Maximum size of an assembly region                          
#    max_reads_per_alignment_start                      Maximum number of reads to retain per alignment start position. Reads above this
#    min_assembly_region_size                           Minimum size of an assembly region                          
#    min_base_quality_score                             Minimum base quality required to consider a base for calling
#    native_pair_hmm_threads                            How many threads should a native pairHMM implementation use 
#    native_pair_hmm_use_double_precision               use double precision in the native pairHmm. This is slower but matches the java 
#    num_reference_samples_if_no_call                   Number of hom-ref genotypes to infer at sites not present in a panel
#    output_mode                                        Specifies which type of calls we should output              
#    pedigree                                           Pedigree file for determining the population "founders"     
#    population_callset                                 Callset to use in calculating genotype priors               
#    sample_name                                        Name of single sample to use from a multi-sample bam        
#    sample_ploidy                                      Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of sa
#    sites_only_vcf_output                              If true, don't emit genotype fields when writing vcf file output.
#    standard_min_confidence_threshold_for_calling      The minimum phred-scaled confidence threshold at which variants should be called
#    version                                            display the version number for this tool                    
#
#  Optional Common Arguments
#    add_output_sam_program_record                      If true, adds a PG tag to created SAM/BAM/CRAM files.       
#    add_output_vcf_command_line                        If true, adds a command line header line to created VCF files.
#    create_output_bam_index                            If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.
#    create_output_bam_md5                              If true, create a MD5 digest for any BAM/SAM/CRAM file created
#    create_output_variant_index                        If true, create a VCF index when writing a coordinate-sorted VCF file.
#    create_output_variant_md5                          If true, create a a MD5 digest any VCF file created.        
#    disable_read_filter                                Read filters to be disabled before analysis                 
#    disable_tool_default_read_filters                  Disable all tool default read filters (WARNING: many tools will not function cor
#    exclude_intervals                                  One or more genomic intervals to exclude from processing    
#    gatk_config_file                                   A configuration file to use with the GATK.                  
#    interval_exclusion_padding                         Amount of padding (in bp) to add to each interval you are excluding.
#    interval_padding                                   Amount of padding (in bp) to add to each interval you are including.
#    interval_set_rule                                  Set merging approach to use for combining interval inputs   
#    lenient                                            Lenient processing of VCF files                             
#    QUIET                                              Whether to suppress job-summary info on System.err.         
#    read_filter                                        Read filters to be applied before analysis                  
#    read_index                                         Indices to use for the read inputs. If specified, an index must be provided for 
#    read_validation_stringency                         Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The 
#    seconds_between_progress_updates                   Output traversal statistics every time this many seconds elapse
#    sequence_dictionary                                Use the given sequence dictionary as the master/canonical sequence dictionary.  
#    tmp_dir                                            Temp directory to use.                                      
#    use_jdk_deflater                                   Whether to use the JdkDeflater (as opposed to IntelDeflater)
#    use_jdk_inflater                                   Whether to use the JdkInflater (as opposed to IntelInflater)
#    verbosity                                          Control verbosity of logging.                               

workflow HaplotypeCaller {

  input {
    #GATK location
    String gatk

    # Required Arguments
    Array[String] I
    String O
    File reference
    File referenceDictionary
    File referenceIndex

    # Optional Tool Arguments
    File? alleles
    Boolean? annotate_with_num_discovered_alleles
    Array[String]? annotation
    Array[String]? annotation_group
    Array[String]? annotations_to_exclude
    Array[File]? arguments_file
    String? assembly_region_out
    Int? assembly_region_padding
    Int? base_quality_score_threshold
    Int? cloud_index_prefetch_buffer
    Int? cloud_prefetch_buffer
    Float? contamination_fraction_to_filter
    Boolean? correct_overlapping_quality
    File? dbsnp
    Boolean? disable_bam_index_caching
    Boolean? disable_sequence_dictionary_validation
    Array[String]? founder_id
    Int? gcs_max_retries
    String? gcs_project_for_requester_pays
    String? graph_output
    Boolean? help
    Float? heterozygosity
    Float? heterozygosity_stdev
    Float? indel_heterozygosity
    String? interval_merging_rule
    Array[String]? intervals
    Int? max_assembly_region_size
    Int? max_reads_per_alignment_start
    Int? min_assembly_region_size
    Int? min_base_quality_score
    Int? native_pair_hmm_threads
    Boolean? native_pair_hmm_use_double_precision
    Int? num_reference_samples_if_no_call
    String? output_mode
    File? pedigree
    File? population_callset
    String? sample_name
    Int? sample_ploidy
    Boolean? sites_only_vcf_output
    Float? standard_min_confidence_threshold_for_calling
    Boolean? version

    # Optional Common Arguments
    Boolean? add_output_sam_program_record
    Boolean? add_output_vcf_command_line
    Boolean? create_output_bam_index
    Boolean? create_output_bam_md5
    Boolean? create_output_variant_index
    Boolean? create_output_variant_md5
    Array[String]? disable_read_filter
    Boolean? disable_tool_default_read_filters
    Array[String]? exclude_intervals
    String? gatk_config_file
    Int? interval_exclusion_padding
    Int? interval_padding
    String? interval_set_rule
    Boolean? lenient
    Boolean? QUIET
    Array[String]? read_filter
    Array[String]? read_index
    String? read_validation_stringency
    Float? seconds_between_progress_updates
    String? sequence_dictionary
    File? tmp_dir
    Boolean? use_jdk_deflater
    Boolean? use_jdk_inflater
    String? verbosity

  }

  call HaplotypeCallerTask {

    input:

    #GATK location
    gatk                                               = gatk,

    # Required Arguments
    I                                                  = I,
    O                                                  = O,
    reference                                          = reference,
    referenceDictionary                                = referenceDictionary,
    referenceIndex                                     = referenceIndex,

    # Optional Tool Arguments
    alleles                                            = alleles,
    annotate_with_num_discovered_alleles               = annotate_with_num_discovered_alleles,
    annotation                                         = annotation,
    annotation_group                                   = annotation_group,
    annotations_to_exclude                             = annotations_to_exclude,
    arguments_file                                     = arguments_file,
    assembly_region_out                                = assembly_region_out,
    assembly_region_padding                            = assembly_region_padding,
    base_quality_score_threshold                       = base_quality_score_threshold,
    cloud_index_prefetch_buffer                        = cloud_index_prefetch_buffer,
    cloud_prefetch_buffer                              = cloud_prefetch_buffer,
    contamination_fraction_to_filter                   = contamination_fraction_to_filter,
    correct_overlapping_quality                        = correct_overlapping_quality,
    dbsnp                                              = dbsnp,
    disable_bam_index_caching                          = disable_bam_index_caching,
    disable_sequence_dictionary_validation             = disable_sequence_dictionary_validation,
    founder_id                                         = founder_id,
    gcs_max_retries                                    = gcs_max_retries,
    gcs_project_for_requester_pays                     = gcs_project_for_requester_pays,
    graph_output                                       = graph_output,
    help                                               = help,
    heterozygosity                                     = heterozygosity,
    heterozygosity_stdev                               = heterozygosity_stdev,
    indel_heterozygosity                               = indel_heterozygosity,
    interval_merging_rule                              = interval_merging_rule,
    intervals                                          = intervals,
    max_assembly_region_size                           = max_assembly_region_size,
    max_reads_per_alignment_start                      = max_reads_per_alignment_start,
    min_assembly_region_size                           = min_assembly_region_size,
    min_base_quality_score                             = min_base_quality_score,
    native_pair_hmm_threads                            = native_pair_hmm_threads,
    native_pair_hmm_use_double_precision               = native_pair_hmm_use_double_precision,
    num_reference_samples_if_no_call                   = num_reference_samples_if_no_call,
    output_mode                                        = output_mode,
    pedigree                                           = pedigree,
    population_callset                                 = population_callset,
    sample_name                                        = sample_name,
    sample_ploidy                                      = sample_ploidy,
    sites_only_vcf_output                              = sites_only_vcf_output,
    standard_min_confidence_threshold_for_calling      = standard_min_confidence_threshold_for_calling,
    version                                            = version,

    # Optional Common Arguments
    add_output_sam_program_record                      = add_output_sam_program_record,
    add_output_vcf_command_line                        = add_output_vcf_command_line,
    create_output_bam_index                            = create_output_bam_index,
    create_output_bam_md5                              = create_output_bam_md5,
    create_output_variant_index                        = create_output_variant_index,
    create_output_variant_md5                          = create_output_variant_md5,
    disable_read_filter                                = disable_read_filter,
    disable_tool_default_read_filters                  = disable_tool_default_read_filters,
    exclude_intervals                                  = exclude_intervals,
    gatk_config_file                                   = gatk_config_file,
    interval_exclusion_padding                         = interval_exclusion_padding,
    interval_padding                                   = interval_padding,
    interval_set_rule                                  = interval_set_rule,
    lenient                                            = lenient,
    QUIET                                              = QUIET,
    read_filter                                        = read_filter,
    read_index                                         = read_index,
    read_validation_stringency                         = read_validation_stringency,
    seconds_between_progress_updates                   = seconds_between_progress_updates,
    sequence_dictionary                                = sequence_dictionary,
    tmp_dir                                            = tmp_dir,
    use_jdk_deflater                                   = use_jdk_deflater,
    use_jdk_inflater                                   = use_jdk_inflater,
    verbosity                                          = verbosity,

  }

  output {
    # Workflow Outputs                                  
    File HaplotypeCallerresults = HaplotypeCallerTask.HaplotypeCallerTask_results
  }
}

task HaplotypeCallerTask {

  input {
    String gatk
    Array[String] I
    String O
    File reference
    File referenceDictionary
    File referenceIndex
    File? alleles
    Boolean? annotate_with_num_discovered_alleles
    Array[String]? annotation
    Array[String]? annotation_group
    Array[String]? annotations_to_exclude
    Array[File]? arguments_file
    String? assembly_region_out
    Int? assembly_region_padding
    Int? base_quality_score_threshold
    Int? cloud_index_prefetch_buffer
    Int? cloud_prefetch_buffer
    Float? contamination_fraction_to_filter
    Boolean? correct_overlapping_quality
    File? dbsnp
    Boolean? disable_bam_index_caching
    Boolean? disable_sequence_dictionary_validation
    Array[String]? founder_id
    Int? gcs_max_retries
    String? gcs_project_for_requester_pays
    String? graph_output
    Boolean? help
    Float? heterozygosity
    Float? heterozygosity_stdev
    Float? indel_heterozygosity
    String? interval_merging_rule
    Array[String]? intervals
    Int? max_assembly_region_size
    Int? max_reads_per_alignment_start
    Int? min_assembly_region_size
    Int? min_base_quality_score
    Int? native_pair_hmm_threads
    Boolean? native_pair_hmm_use_double_precision
    Int? num_reference_samples_if_no_call
    String? output_mode
    File? pedigree
    File? population_callset
    String? sample_name
    Int? sample_ploidy
    Boolean? sites_only_vcf_output
    Float? standard_min_confidence_threshold_for_calling
    Boolean? version
    Boolean? add_output_sam_program_record
    Boolean? add_output_vcf_command_line
    Boolean? create_output_bam_index
    Boolean? create_output_bam_md5
    Boolean? create_output_variant_index
    Boolean? create_output_variant_md5
    Array[String]? disable_read_filter
    Boolean? disable_tool_default_read_filters
    Array[String]? exclude_intervals
    String? gatk_config_file
    Int? interval_exclusion_padding
    Int? interval_padding
    String? interval_set_rule
    Boolean? lenient
    Boolean? QUIET
    Array[String]? read_filter
    Array[String]? read_index
    String? read_validation_stringency
    Float? seconds_between_progress_updates
    String? sequence_dictionary
    File? tmp_dir
    Boolean? use_jdk_deflater
    Boolean? use_jdk_inflater
    String? verbosity

  }

  command <<<
    ~{gatk} HaplotypeCaller \
    --I ~{sep=' --I ' I} \
    --O ~{sep=' --O ' O} \
    --reference ~{sep=' --reference ' reference} \
    ~{true='--alleles ' false='' defined(alleles)}~{sep=' --alleles ' alleles} \
    ~{true='--annotate_with_num_discovered_alleles ' false='' defined(annotate_with_num_discovered_alleles)}~{sep=' --annotate_with_num_discovered_alleles ' annotate_with_num_discovered_alleles} \
    ~{true='--annotation ' false='' defined(annotation)}~{sep=' --annotation ' annotation} \
    ~{true='--annotation_group ' false='' defined(annotation_group)}~{sep=' --annotation_group ' annotation_group} \
    ~{true='--annotations_to_exclude ' false='' defined(annotations_to_exclude)}~{sep=' --annotations_to_exclude ' annotations_to_exclude} \
    ~{true='--arguments_file ' false='' defined(arguments_file)}~{sep=' --arguments_file ' arguments_file} \
    ~{true='--assembly_region_out ' false='' defined(assembly_region_out)}~{sep=' --assembly_region_out ' assembly_region_out} \
    ~{true='--assembly_region_padding ' false='' defined(assembly_region_padding)}~{sep=' --assembly_region_padding ' assembly_region_padding} \
    ~{true='--base_quality_score_threshold ' false='' defined(base_quality_score_threshold)}~{sep=' --base_quality_score_threshold ' base_quality_score_threshold} \
    ~{true='--cloud_index_prefetch_buffer ' false='' defined(cloud_index_prefetch_buffer)}~{sep=' --cloud_index_prefetch_buffer ' cloud_index_prefetch_buffer} \
    ~{true='--cloud_prefetch_buffer ' false='' defined(cloud_prefetch_buffer)}~{sep=' --cloud_prefetch_buffer ' cloud_prefetch_buffer} \
    ~{true='--contamination_fraction_to_filter ' false='' defined(contamination_fraction_to_filter)}~{sep=' --contamination_fraction_to_filter ' contamination_fraction_to_filter} \
    ~{true='--correct_overlapping_quality ' false='' defined(correct_overlapping_quality)}~{sep=' --correct_overlapping_quality ' correct_overlapping_quality} \
    ~{true='--dbsnp ' false='' defined(dbsnp)}~{sep=' --dbsnp ' dbsnp} \
    ~{true='--disable_bam_index_caching ' false='' defined(disable_bam_index_caching)}~{sep=' --disable_bam_index_caching ' disable_bam_index_caching} \
    ~{true='--disable_sequence_dictionary_validation ' false='' defined(disable_sequence_dictionary_validation)}~{sep=' --disable_sequence_dictionary_validation ' disable_sequence_dictionary_validation} \
    ~{true='--founder_id ' false='' defined(founder_id)}~{sep=' --founder_id ' founder_id} \
    ~{true='--gcs_max_retries ' false='' defined(gcs_max_retries)}~{sep=' --gcs_max_retries ' gcs_max_retries} \
    ~{true='--gcs_project_for_requester_pays ' false='' defined(gcs_project_for_requester_pays)}~{sep=' --gcs_project_for_requester_pays ' gcs_project_for_requester_pays} \
    ~{true='--graph_output ' false='' defined(graph_output)}~{sep=' --graph_output ' graph_output} \
    ~{true='--help ' false='' defined(help)}~{sep=' --help ' help} \
    ~{true='--heterozygosity ' false='' defined(heterozygosity)}~{sep=' --heterozygosity ' heterozygosity} \
    ~{true='--heterozygosity_stdev ' false='' defined(heterozygosity_stdev)}~{sep=' --heterozygosity_stdev ' heterozygosity_stdev} \
    ~{true='--indel_heterozygosity ' false='' defined(indel_heterozygosity)}~{sep=' --indel_heterozygosity ' indel_heterozygosity} \
    ~{true='--interval_merging_rule ' false='' defined(interval_merging_rule)}~{sep=' --interval_merging_rule ' interval_merging_rule} \
    ~{true='--intervals ' false='' defined(intervals)}~{sep=' --intervals ' intervals} \
    ~{true='--max_assembly_region_size ' false='' defined(max_assembly_region_size)}~{sep=' --max_assembly_region_size ' max_assembly_region_size} \
    ~{true='--max_reads_per_alignment_start ' false='' defined(max_reads_per_alignment_start)}~{sep=' --max_reads_per_alignment_start ' max_reads_per_alignment_start} \
    ~{true='--min_assembly_region_size ' false='' defined(min_assembly_region_size)}~{sep=' --min_assembly_region_size ' min_assembly_region_size} \
    ~{true='--min_base_quality_score ' false='' defined(min_base_quality_score)}~{sep=' --min_base_quality_score ' min_base_quality_score} \
    ~{true='--native_pair_hmm_threads ' false='' defined(native_pair_hmm_threads)}~{sep=' --native_pair_hmm_threads ' native_pair_hmm_threads} \
    ~{true='--native_pair_hmm_use_double_precision ' false='' defined(native_pair_hmm_use_double_precision)}~{sep=' --native_pair_hmm_use_double_precision ' native_pair_hmm_use_double_precision} \
    ~{true='--num_reference_samples_if_no_call ' false='' defined(num_reference_samples_if_no_call)}~{sep=' --num_reference_samples_if_no_call ' num_reference_samples_if_no_call} \
    ~{true='--output_mode ' false='' defined(output_mode)}~{sep=' --output_mode ' output_mode} \
    ~{true='--pedigree ' false='' defined(pedigree)}~{sep=' --pedigree ' pedigree} \
    ~{true='--population_callset ' false='' defined(population_callset)}~{sep=' --population_callset ' population_callset} \
    ~{true='--sample_name ' false='' defined(sample_name)}~{sep=' --sample_name ' sample_name} \
    ~{true='--sample_ploidy ' false='' defined(sample_ploidy)}~{sep=' --sample_ploidy ' sample_ploidy} \
    ~{true='--sites_only_vcf_output ' false='' defined(sites_only_vcf_output)}~{sep=' --sites_only_vcf_output ' sites_only_vcf_output} \
    ~{true='--standard_min_confidence_threshold_for_calling ' false='' defined(standard_min_confidence_threshold_for_calling)}~{sep=' --standard_min_confidence_threshold_for_calling ' standard_min_confidence_threshold_for_calling} \
    ~{true='--version ' false='' defined(version)}~{sep=' --version ' version} \
    ~{true='--add_output_sam_program_record ' false='' defined(add_output_sam_program_record)}~{sep=' --add_output_sam_program_record ' add_output_sam_program_record} \
    ~{true='--add_output_vcf_command_line ' false='' defined(add_output_vcf_command_line)}~{sep=' --add_output_vcf_command_line ' add_output_vcf_command_line} \
    ~{true='--create_output_bam_index ' false='' defined(create_output_bam_index)}~{sep=' --create_output_bam_index ' create_output_bam_index} \
    ~{true='--create_output_bam_md5 ' false='' defined(create_output_bam_md5)}~{sep=' --create_output_bam_md5 ' create_output_bam_md5} \
    ~{true='--create_output_variant_index ' false='' defined(create_output_variant_index)}~{sep=' --create_output_variant_index ' create_output_variant_index} \
    ~{true='--create_output_variant_md5 ' false='' defined(create_output_variant_md5)}~{sep=' --create_output_variant_md5 ' create_output_variant_md5} \
    ~{true='--disable_read_filter ' false='' defined(disable_read_filter)}~{sep=' --disable_read_filter ' disable_read_filter} \
    ~{true='--disable_tool_default_read_filters ' false='' defined(disable_tool_default_read_filters)}~{sep=' --disable_tool_default_read_filters ' disable_tool_default_read_filters} \
    ~{true='--exclude_intervals ' false='' defined(exclude_intervals)}~{sep=' --exclude_intervals ' exclude_intervals} \
    ~{true='--gatk_config_file ' false='' defined(gatk_config_file)}~{sep=' --gatk_config_file ' gatk_config_file} \
    ~{true='--interval_exclusion_padding ' false='' defined(interval_exclusion_padding)}~{sep=' --interval_exclusion_padding ' interval_exclusion_padding} \
    ~{true='--interval_padding ' false='' defined(interval_padding)}~{sep=' --interval_padding ' interval_padding} \
    ~{true='--interval_set_rule ' false='' defined(interval_set_rule)}~{sep=' --interval_set_rule ' interval_set_rule} \
    ~{true='--lenient ' false='' defined(lenient)}~{sep=' --lenient ' lenient} \
    ~{true='--QUIET ' false='' defined(QUIET)}~{sep=' --QUIET ' QUIET} \
    ~{true='--read_filter ' false='' defined(read_filter)}~{sep=' --read_filter ' read_filter} \
    ~{true='--read_index ' false='' defined(read_index)}~{sep=' --read_index ' read_index} \
    ~{true='--read_validation_stringency ' false='' defined(read_validation_stringency)}~{sep=' --read_validation_stringency ' read_validation_stringency} \
    ~{true='--seconds_between_progress_updates ' false='' defined(seconds_between_progress_updates)}~{sep=' --seconds_between_progress_updates ' seconds_between_progress_updates} \
    ~{true='--sequence_dictionary ' false='' defined(sequence_dictionary)}~{sep=' --sequence_dictionary ' sequence_dictionary} \
    ~{true='--tmp_dir ' false='' defined(tmp_dir)}~{sep=' --tmp_dir ' tmp_dir} \
    ~{true='--use_jdk_deflater ' false='' defined(use_jdk_deflater)}~{sep=' --use_jdk_deflater ' use_jdk_deflater} \
    ~{true='--use_jdk_inflater ' false='' defined(use_jdk_inflater)}~{sep=' --use_jdk_inflater ' use_jdk_inflater} \
    ~{true='--verbosity ' false='' defined(verbosity)}~{sep=' --verbosity ' verbosity} \
  >>>


  output {
    # Task Outputs                                      
    File HaplotypeCallerTask_results = stdout()
  }
 }

