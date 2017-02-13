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
     */
    public final static String DOC_CAT_SPARK = "Spark Tools";
    public final static String DOC_CAT_SPARK_SUMMARY = "Tools that use Apache Spark for scaling out (experimental)";

    public final static String DOC_CAT_SPARK_PIPELINE = "Spark Pipelines";
    public final static String DOC_CAT_SPARK_PIPELINE_SUMMARY = "Pipelines that combine tools and use Apache Spark for scaling out (experimental)";

    public final static String DOC_CAT_READFILTERS = "Read Filters";
    public final static String DOC_CAT_READFILTERS_SUMMARY = "Read Filters used by the engine to select reads to be included for analysis";

    public final static String DOC_CAT_FASTA = "Fasta File Tools";
    public final static String DOC_CAT_FASTA_SUMMARY = "Tools for analysis and manipulation of files in fasta format";

    public final static String DOC_CAT_INTERVALS = "Interval Tools";
    public final static String DOC_CAT_INTERVALS_SUMMARY = "Tools for processing intervals and associated overlapping records";

    public final static String DOC_CAT_READS = "SAM/BAM/CRAM Tools";
    public final static String DOC_CAT_READS_SUMMARY = "Tools for manipulating read-level data (SAM/BAM/CRAM)";

    public final static String DOC_CAT_QC = "Diagnostics and Quality Control Tools";
    public final static String DOC_CAT_QC_SUMMARY = "Tools for Diagnostics and Quality Control";

    public final static String DOC_CAT_CNV = "Copy Number Analysis Tools";
    public final static String DOC_CAT_CNV_SUMMARY = "Tools to analyze copy number data.";

    public final static String DOC_CAT_SPARK_SV = "Spark Structural Variation Analysis Tools";
    public final static String DOC_CAT_SPARK_SV_SUMMARY = "Structural variation analysis tools that use Apache Spark for scaling out (experimental)";

    public final static String DOC_CAT_VARIANT = "VCF Tools";
    public final static String DOC_CAT_VARIANT_SUMMARY = "Tools for manipulating variants and associated metadata";

    public final static String DOC_CAT_SPARK_VALIDATION = "Spark Test Tools";
    public final static String DOC_CAT_SPARK_VALIDATION_SUMMARY = "Tools for experimenting with Spark";

    public final static String DOC_CAT_TEST = "Test Tools";
    public final static String DOC_CAT_TEST_SUMMARY = "Tools for internal test purposes";

    public final static String DOC_CAT_VARDISC = "Variant Discovery Tools";
    public final static String DOC_CAT_VAREVAL = "Variant Evaluation Tools";
    public final static String DOC_CAT_VARMANIP = "Variant Manipulation Tools";

    public final static String DOC_CAT_EXAMPLE = "Examples (Exclude)";

    /**
     * List of "supercategory" values used for doc purposes. Every doc group name can/should be put into
     * one of the following supercategories.
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
            groupToSuperCategory.put(DOC_CAT_SPARK, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SPARK_PIPELINE, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SPARK_SV, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_FASTA, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_READS, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_QC, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_CNV, DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_VARIANT, DOC_SUPERCAT_TOOLS);

            // supercat Utilities
            groupToSuperCategory.put(DOC_CAT_READFILTERS, DOC_SUPERCAT_UTILITIES);

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