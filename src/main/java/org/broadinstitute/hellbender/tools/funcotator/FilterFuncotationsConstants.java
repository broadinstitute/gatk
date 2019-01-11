package org.broadinstitute.hellbender.tools.funcotator;

public class FilterFuncotationsConstants {
    /**
     * Key for the INFO field added to all variants by {@link FilterFuncotations},
     * indicating the clinical significance (if any) of the Funcotations on that variant.
     */
    public static final String CLINSIG_INFO_KEY = "CLINSIG";

    /**
     * Description for {@value CLINSIG_INFO_KEY} to include in VCF headers.
     */
    public static final String CLINSIG_INFO_KEY_DESCRIPTION =
            "Rule(s) which caused this annotation to be flagged as clinically significant.";

    /**
     * Value to assign to {@value CLINSIG_INFO_KEY} for variants that have no
     * clinically-significant Funcotations.
     */
    public static final String CLINSIG_INFO_NOT_SIGNIFICANT = "NONE";

    /**
     * FILTER value applied by {@link FilterFuncotations} to all variants which have
     * no clinically-significant Funcotations.
     */
    public static final String NOT_CLINSIG_FILTER = "NOT_" + CLINSIG_INFO_KEY;

    /**
     * Description for {@value NOT_CLINSIG_FILTER} to include in VCF headers.
     */
    public static final String NOT_CLINSIG_FILTER_DESCRIPTION = "Filter for clinically insignificant variants.";

    /**
     * Delimiting string to place between values in the {@value CLINSIG_INFO_KEY} INFO field
     * when Funcotations for a variant match multiple filters.
     */
    public static final String FILTER_DELIMITER = ",";
}
