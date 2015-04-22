package org.broadinstitute.hellbender.utils.help;

public final class HelpConstants {

    public final static String BASE_GATK_URL = "http://www.broadinstitute.org/gatk";
    public final static String GATK_DOCS_URL = BASE_GATK_URL + "/tooldocs/";
    public final static String GATK_FORUM_URL = "http://gatkforums.broadinstitute.org/";
    public final static String GATK_FORUM_API_URL = "https://gatkforums.broadinstitute.org/api/v1/";

    /**
     * Definition of the group names / categories of tools.
     * The names get parsed to make supercategories in the doc index,
     * so be careful when making big changes -- see GATKDoclet.java toMap()
     */
    public final static String DOCS_CAT_DATA = "Sequence Data Processing Tools";
    public final static String DOCS_CAT_QC = "Diagnostics and Quality Control Tools";
    public final static String DOCS_CAT_ENGINE = "Engine Parameters (available to all tools)";
    public final static String DOCS_CAT_RF = "Read Filters";
    public final static String DOCS_CAT_REFUTILS = "Reference Utilities";
    public final static String DOCS_CAT_RODCODECS = "ROD Codecs";
    public final static String DOCS_CAT_USRERR = "User Exceptions (DevZone)";
    public final static String DOCS_CAT_VALIDATION = "Validation Utilities";
    public final static String DOCS_CAT_ANNOT = "Variant Annotations";
    public final static String DOCS_CAT_VARDISC = "Variant Discovery Tools";
    public final static String DOCS_CAT_VARMANIP = "Variant Evaluation and Manipulation Tools";
    public final static String DOCS_CAT_TOY = "Toy Walkers (DevZone)";
    public final static String DOCS_CAT_HELPUTILS = "Help Utilities";

    public static String forumPost(String post) {
        return GATK_FORUM_URL + post;
    }

}