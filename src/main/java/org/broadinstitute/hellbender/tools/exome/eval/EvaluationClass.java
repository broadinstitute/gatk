package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Enumeration of possible evaluation case outcomes or classes such as {@link #TRUE_POSITIVE}s,
 * {@link #FALSE_NEGATIVE} and so forth.
 * <p>
 * Includes some utility methods and fields to deal with VCF representation of these evaluation classes.
 * </p>
 */
enum EvaluationClass {

    /**
     * When truth and calls overlap and are all compatible (all are deletion or all are duplication).
     */
    TRUE_POSITIVE("TP", "When truth and calls overlap and are all compatible (all are deletion or all are duplication)"),

    /**
     * When truth does not overlap with any call.
     */
    FALSE_NEGATIVE("FN", "When truth does not overlap with any call"),

    /**
     * When truth overlaps with several calls that are discordant amongst them.
     */
    MIXED_POSITIVE("MP", "When truth overlaps with several calls that are discordant amongst them"),

    /**
     * When truth overlaps with several calls that are concordant amongst them but
     * the are discordant with the truth.
     */
    DISCORDANT_POSITIVE("DP", "When truth overlaps with several calls that are concordant amongst them but the are discordant with the truth"),

    /**
     * When a call does not overlap with any truth call.
     */
    UNKNOWN_POSITIVE("UP", "When a call does not overlap with any truth call"),

    /**
     * When a call overlaps with a ref truth call.
     */
    FALSE_POSITIVE("FP", "When a call overlap with a ref truth call" );

    /**
     * The key for vcf header lines describing evaluation classes.
     */
    public static final String VCF_HEADER_KEY = "EVAL_CLASS";

    private static final Map<String, EvaluationClass> VALUE_BY_ACRONYM =
            Collections.unmodifiableMap(Stream.of(values())
                    .collect(Collectors.toMap(s -> s.acronym, s -> s)));

    /**
     * Short acronym name for the class.
     */
    public final String acronym;

    /**
     * Description of this class meaning.
     */
    public final String description;

    /**
     * Creates a new evaluation class type given its acronym.
     * @param acronym the evaluation class acronym text.
     * @param description the evaluation class description text.
     * @throws IllegalArgumentException if either {@code acronym} or {@code description} is {@code null}.
     */
    EvaluationClass(final String acronym, final String description) {
        this.acronym = Utils.nonNull(acronym);
        this.description = Utils.nonNull(description);
    }

    /**
     * Represents the evaluation class into a string using its acronym.
     * @return never {@code null}.
     */
    @Override
    public String toString() {
        return acronym;
    }

    /**
     * Composes the vcf-header line to describe this evaluation class using its acronym as the ID.
     * @return never {@code null}.
     */
    public VCFHeaderLine headerLine() {
        return new VCFSimpleHeaderLine(VCF_HEADER_KEY, acronym, description);
    }

    /**
     * Add all pertinent meta-data header lines to a vcf header.
     * @param header the target header.
     * @throws IllegalArgumentException if {@code header} is {@code null}.
     */
    public static void addHeaderLinesTo(final VCFHeader header) {
        for (final EvaluationClass value : values()) {
            header.addMetaDataLine(value.headerLine());
        }
    }

    /**
     * Resolves the evaluation class given either is acronym or full name.
     * <p>
     * This method should be used in place of {@link #valueOf} when one wants to
     * accept acronym as a valid string representation of an evaluation class.
     * </p>
     *
     * @param string the query string to parse into a evaluation class.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code string} is {@code null} or there is no such a evaluation class
     */
    public static EvaluationClass parseString(final String string) {
        if (VALUE_BY_ACRONYM.containsKey(string)) {
            return VALUE_BY_ACRONYM.get(string);
        } else {
            return valueOf(string);
        }
    }
}
