package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public final class HelpConstants {

    private HelpConstants() {};

    public static final String GATK_FORUM_URL = "http://gatkforums.broadinstitute.org/";

    public static String forumPost(String post) {
        return GATK_FORUM_URL + post;
    }

    /**
     * Definition of the group names / descriptions for documentation/help purposes.
     *
     */

    // Start temporary Picard placeholders. The Program Groups that correspond to these strings should be
    // replaced with references to Picard classes when they are available.
    public final static String DOC_CAT_BAM_PREPROCESSING = "Alignment, duplicate flagging and BQSR"; // Picard placeholder
    public final static String DOC_CAT_BAM_PREPROCESSING_SUMMARY = "Tools that align reads, flag duplicates and recalibrate base qualities";

    public final static String DOC_CAT_DIAGNOSTICS_AND_QC = "Diagnostics and Quality Control"; // Picard placeholder
    public final static String DOC_CAT_DIAGNOSTICS_AND_QC_SUMMARY = "Tools that collect sequencing quality-related and comparative metrics";

    public final static String DOC_CAT_INTERVALS = "Intervals Manipulation"; // Picard placeholder
    public final static String DOC_CAT_INTERVALS_SUMMARY = "Tools that process genomic intervals in various formats";

    public final static String DOC_CAT_OTHER = "Other"; // Picard placeholder
    public final static String DOC_CAT_OTHER_SUMMARY = "Miscellaneous tools, e.g. those that aid in data streaming";

    public final static String DOC_CAT_READ_DATA = "Read Data"; // Picard placeholder
    public final static String DOC_CAT_READ_DATA_SUMMARY = "Tools that manipulate read data in SAM, BAM or CRAM format";

    public final static String DOC_CAT_REFERENCE = "Reference"; // Picard placeholder
    public final static String DOC_CAT_REFERENCE_SUMMARY = "Tools that analyze and manipulate FASTA format references";

    public final static String DOC_CAT_VARIANT = "VCF Tools";  // Picard placeholder
    public final static String DOC_CAT_VARIANT_SUMMARY = "Tools for manipulating variants and associated metadata";

    public final static String DOC_CAT_VARIANT_FILTERING = "Variant Filtering"; // Picard placeholder
    public final static String DOC_CAT_VARIANT_FILTERING_SUMMARY = "Tools that filter variants";

    public final static String DOC_CAT_VARIANT_EVALUATION = "Variant Evaluation and Refinement";   // Picard placeholder
    public final static String DOC_CAT_VARIANT_EVALUATION_SUMMARY = "Tools that evaluate and refine variant calls, e.g. by adding annotations that are not calculated during variant calling";

    public final static String DOC_CAT_VCF_MANIPULATION = "VCF Manipulation";   // Picard placeholder
    public final static String DOC_CAT_VCF_MANIPULATION_SUMMARY = "Tools that manipulate variant call format (VCF) data";
    // End temporary Picard placeholders.

    // Start GATK Program groups
    public final static String DOC_CAT_ANNOTATORS = "Variant Annotations";
    public final static String DOC_CAT_ANNOTATORS_SUMMARY = "Available to HaplotypeCaller, Mutect2, VariantAnnotator and GenotypeGVCFs. " +
            "See https://software.broadinstitute.org/gatk/documentation/article?id=10836";

    public final static String DOC_CAT_COVERAGE_ANALYSIS = "Coverage Analysis";
    public final static String DOC_CAT_COVERAGE_ANALYSIS_SUMMARY = "Tools that count coverage, e.g. depth per allele";

    public final static String DOC_CAT_CNV = "Copy Number Variant Discovery";
    public final static String DOC_CAT_CNV_SUMMARY = "Tools that analyze read coverage to detect copy number variants.";

    public final static String DOC_CAT_EXAMPLE = "Example Tools";
    public final static String DOC_CAT_EXAMPLE_SUMMARY = "Example tools that show developers how to implement new tools";

    public final static String DOC_CAT_METAGENOMICS = "Metagenomics";
    public final static String DOC_CAT_METAGENOMICS_SUMMARY = "Tools that perform metagenomic analysis, e.g. microbial community composition and pathogen detection";

    public final static String DOC_CAT_QC = "Diagnostics and Quality Control Tools";
    public final static String DOC_CAT_QC_SUMMARY = "Tools for Diagnostics and Quality Control";

    public final static String DOC_CAT_READFILTERS = "Read Filters";
    public final static String DOC_CAT_READFILTERS_SUMMARY = "Applied by engine to select reads for analysis";

    public final static String DOC_CAT_SHORT_VARIANT_DISCOVERY = "Short Variant Discovery";
    public final static String DOC_CAT_SHORT_VARIANT_DISCOVERY_SUMMARY = "Tools that perform variant calling and genotyping for short variants (SNPs, SNVs and Indels)";

    public final static String DOC_CAT_SV_DISCOVERY = "Structural Variant Discovery";
    public final static String DOC_CAT_SV_DISCOVERY_SUMMARY = "Tools that detect structural variants";

    public final static String DOC_CAT_SPARK_SV = "Spark Structural Variation Analysis Tools";
    public final static String DOC_CAT_SPARK_SV_SUMMARY = "Structural variation analysis tools that use Apache Spark for scaling out (experimental)";
    // End GATK Program groups

    // Obsolete - to be removed
    public final static String DOC_CAT_BWA_MEM_UTILS = "Bwa Mem JNI binding Tools";
    public final static String DOC_CAT_BWA_MEM_UTILS_SUMMARY = "Bwa Mem JNI binding Tools";

    public final static String DOC_CAT_SPARK_PATHSEQ = "PathSeq Tools";
    public final static String DOC_CAT_SPARK_PATHSEQ_SUMMARY = "Tools for detecting pathogens (experimental)";

    public final static String DOC_CAT_SPARK_VALIDATION = "Spark Test Tools";
    public final static String DOC_CAT_SPARK_VALIDATION_SUMMARY = "Tools for experimenting with Spark";

    public final static String DOC_CAT_SPARK = "Spark Tools";
    public final static String DOC_CAT_SPARK_SUMMARY = "Tools that use Apache Spark for scaling out (experimental)";

    public final static String DOC_CAT_SPARK_PIPELINE = "Spark Pipelines";
    public final static String DOC_CAT_SPARK_PIPELINE_SUMMARY = "Pipelines that combine tools and use Apache Spark for scaling out (experimental)";

    public final static String DOC_CAT_TEST = "Test Tools";
    public final static String DOC_CAT_TEST_SUMMARY = "Tools for internal test purposes";

    public final static String DOC_CAT_RNA = "RNA-Specific Tools";
    public final static String DOC_CAT_RNA_SUMMARY = "Tools intended to be used for processing RNA data.";

    /**
     * List of "supercategory" values used for doc purposes. Every doc group name can/should be put into
     * one of the following super-categories.
     */
    private final static String DOC_SUPERCAT_TOOLS = "tools";
    private final static String DOC_SUPERCAT_UTILITIES = "utilities";
    private final static String DOC_SUPERCAT_EXCLUDE = "exclude";

    private static Map<String, String> groupToSuperCategory;
    private static Set<String> sparkCategory;

    private static Map<String, String> getSuperCategoryMap() {
        if (groupToSuperCategory == null) {

            // do this only on demand since we only need it during docgen
            groupToSuperCategory = new HashMap<>();

            // supercat Tools

            // start Picard placeholders
            groupToSuperCategory.put(DOC_CAT_DIAGNOSTICS_AND_QC, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_INTERVALS, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_OTHER, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_READ_DATA, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_REFERENCE, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_READ_DATA, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_VARIANT, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_VARIANT_FILTERING, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_VARIANT_EVALUATION, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_VCF_MANIPULATION, DOC_SUPERCAT_TOOLS);
            // end Picard placeholders

            groupToSuperCategory.put(DOC_CAT_QC, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_CNV, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_COVERAGE_ANALYSIS, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_METAGENOMICS, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SHORT_VARIANT_DISCOVERY, DOC_SUPERCAT_TOOLS);

            groupToSuperCategory.put(DOC_CAT_SPARK, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SPARK_PIPELINE, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SPARK_SV, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SV_DISCOVERY, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SPARK_PATHSEQ, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SPARK_VALIDATION, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_TEST, DOC_SUPERCAT_TOOLS);

            // supercat Utilities
            groupToSuperCategory.put(DOC_CAT_READFILTERS, DOC_SUPERCAT_UTILITIES);
            groupToSuperCategory.put(DOC_CAT_ANNOTATORS, DOC_SUPERCAT_UTILITIES);

            // supercat Exclude
            groupToSuperCategory.put(DOC_CAT_EXAMPLE, DOC_SUPERCAT_EXCLUDE);
            groupToSuperCategory.put(DOC_CAT_SPARK_VALIDATION, DOC_SUPERCAT_EXCLUDE);
            groupToSuperCategory.put(DOC_CAT_TEST, DOC_SUPERCAT_EXCLUDE);
        }
        return groupToSuperCategory;
    }

    /**
     * Given a group name, return a supercategory string for use by the online doc system to determine which
     * supercateogry the group is in. The strings returned by this method should match those used in the
     * corresponding help template.
     *
     * @param groupName
     * @return supercategory string corresponding to {@code groupName} for use in for use in determining
     * which supercateogry the group is in for online doc purposes.
     */
    public static String getSuperCategoryProperty(final String groupName) {
        Utils.nonNull(groupName);
        return getSuperCategoryMap().getOrDefault(groupName, "other");
    }

    private static Set<String> getSparkCategoryMap() {
        // do this on demand since we only need it during docgen
        if (sparkCategory == null) {
            sparkCategory = new HashSet<>();

            sparkCategory.add(DOC_CAT_SPARK);
            sparkCategory.add(DOC_CAT_SPARK_PIPELINE);
            sparkCategory.add(DOC_CAT_SPARK_SV);
        }
        return sparkCategory;
    }

    /**
     * Given a group name, return a string for use by the online doc system to determine if the group is
     * categorized as a "spark" group. The strings returned by this method should match those used in the
     * corresponding help template.
     *
     * @param groupName
     * @return spark category string corresponding to {@code groupName} for use in online doc
     */
    public static String getSparkCategoryProperty(final String groupName) {
        Utils.nonNull(groupName);
        return getSparkCategoryMap().contains(groupName) ? "spark" : "other";
    }

}
