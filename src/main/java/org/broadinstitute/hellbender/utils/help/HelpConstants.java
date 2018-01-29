package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public final class HelpConstants {

    private HelpConstants() {};

    public static final String GATK_FORUM_URL = "http://gatkforums.broadinstitute.org/";
    public static final String GATK_MAIN_SITE = "https://software.broadinstitute.org/gatk/";
    
    public static String forumPost(String post) {
        return GATK_FORUM_URL + post;
    }

    /**
     * Definition of the group names / descriptions for documentation/help purposes.
     *
     */

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

    public final static String DOC_CAT_READFILTERS = "Read Filters";
    public final static String DOC_CAT_READFILTERS_SUMMARY = "Applied by engine to select reads for analysis";

    public final static String DOC_CAT_SHORT_VARIANT_DISCOVERY = "Short Variant Discovery";
    public final static String DOC_CAT_SHORT_VARIANT_DISCOVERY_SUMMARY = "Tools that perform variant calling and genotyping for short variants (SNPs, SNVs and Indels)";

    public final static String DOC_CAT_SV_DISCOVERY = "Structural Variant Discovery";
    public final static String DOC_CAT_SV_DISCOVERY_SUMMARY = "Tools that detect structural variants";

    // End GATK Program groups

    private static Map<String, String> groupToSuperCategory;

    public static Map<String, String> getSuperCategoryMap() {
        if (groupToSuperCategory == null) {

            // do this only on demand since we only need it during docgen

            // start with Picard program groups
            groupToSuperCategory = picard.util.help.HelpConstants.getSuperCategoryMap();

            // supercat Tools
            groupToSuperCategory.put(DOC_CAT_CNV, picard.util.help.HelpConstants.DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_COVERAGE_ANALYSIS, picard.util.help.HelpConstants.DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_METAGENOMICS, picard.util.help.HelpConstants.DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SHORT_VARIANT_DISCOVERY, picard.util.help.HelpConstants.DOC_SUPERCAT_TOOLS);
            groupToSuperCategory.put(DOC_CAT_SV_DISCOVERY, picard.util.help.HelpConstants.DOC_SUPERCAT_TOOLS);

            // supercat Utilities
            groupToSuperCategory.put(DOC_CAT_READFILTERS, picard.util.help.HelpConstants.DOC_SUPERCAT_UTILITIES);
            groupToSuperCategory.put(DOC_CAT_ANNOTATORS, picard.util.help.HelpConstants.DOC_SUPERCAT_UTILITIES);

            // supercat Exclude
            groupToSuperCategory.put(DOC_CAT_EXAMPLE, picard.util.help.HelpConstants.DOC_SUPERCAT_EXCLUDE);
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

}
