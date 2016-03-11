package org.broadinstitute.hellbender.utils.hmm;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Three-state copy number hidden states enum.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public enum CopyNumberTriState {

    /**
     * Represents a loss of one or more copies in a region.
     */
    DELETION("-", Allele.create("<DEL>", false)),

    /**
     * Represents not a loss nor a gain of any copy in a region.
     */
    NEUTRAL("0", Allele.create("N", true)),

    /**
     * Represents a gain of one or more copies in a region.
     */
    DUPLICATION("+", Allele.create("<DUP>", false));

    /**
     * Unmodifiable list alleles where the first allele is the reference allele
     * (corresponding {@link #NEUTRAL}). Other states alternative alleles follow
     * in ascending order based on their ordinal.
     */
    public static List<Allele> ALL_ALLELES =
            Collections.unmodifiableList(Arrays.asList(NEUTRAL.allele, DELETION.allele, DUPLICATION.allele));

    /**
     * Unmodifiable list of alternative alleles order by their corresponding state ordinal.
     */
    public static List<Allele> ALTERNATIVE_ALLELES =
            ALL_ALLELES.subList(1, ALL_ALLELES.size());

    public final String callString;

    public final Allele allele;

    CopyNumberTriState(final String callString, final Allele allele) {
        this.callString = callString;
        this.allele = allele;
    }

    /**
     * The string representation of the call.
     *
     * @return never {@code null}.
     */
    public String toCallString() {
        return callString;
    }

    /**
     * Returns the state given a call string.
     * @param callString the query call string.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code callString} is either a {@code null} or an invalid
     *  call string.
     */
    public static CopyNumberTriState fromCallString(final String callString) {
        Utils.nonNull(callString);
        if (callString.equals(DELETION.callString)) {
            return DELETION;
        } else if (callString.equals(NEUTRAL.callString)) {
            return NEUTRAL;
        } else if (callString.equals(DUPLICATION.callString)) {
            return DUPLICATION;
        } else {
            throw new IllegalArgumentException(String.format("'%s' is not a valid call name", callString));
        }
    }
}
