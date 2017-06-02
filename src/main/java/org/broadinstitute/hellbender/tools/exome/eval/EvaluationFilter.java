package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.stream.Stream;

/**
 * Filters applied to individual truth or called segments when taking them into account for evaluation.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
enum EvaluationFilter {

    /**
     * The segment has a low quality value.
     */
    LowQuality("LQ", "Low Quality"),

    /**
     * The segment is to short to consider.
     */
    ShortEvent("SE", "Short Event"),

    /**
     * The event is common amongst the truth calls.
     */
    CommonEvent("CE", "Common Event"),

    /**
     * The truth is multi-allelic (has deletions and duplications) at that site.
     */
    MultiAllelicTruth("TruthMA", "Multi-allelic CNV Truth"),

    /**
     * The calls are multi-allelic (have deletions and duplications) at that site.
     */
    MultiAllelicCalls("CallsMA", "Multi-allelic CNV Calls"),

    /**
     * There are no qualifying calls at that site.
     */
    NoQualifyingCalls("NQC", "No filter passing calls overlap the truth");

    /**
     * Represent an empty filter set, i.e. the segment pass all filters.
     */
    public static final String PASS = VCFConstants.PASSES_FILTERS_v4;

    /**
     * Shorter name of the filter, used to represent the filter in genotypes.
     */
    public final String acronym;

    /**
     * Short description of the meaning of the filter.
     */
    public final String description;

    EvaluationFilter(final String acronym, final String description) {
        this.acronym = Utils.nonNull(acronym);
        this.description = Utils.nonNull(description);
    }

    /**
     * Convert a string into a filter by matching it to a filter name or acronym.
     *
     * @param str the string to transform.
     * @return {@code null} iff there is such a filter.
     * @throws IllegalArgumentException if {@code str} is {@code null}.
     */
    public static EvaluationFilter fromString(final String str) {
        return Stream.of(values()).filter(f -> f.name().equals(str) || f.acronym.equals(str)).findFirst().orElse(null);
    }
}
