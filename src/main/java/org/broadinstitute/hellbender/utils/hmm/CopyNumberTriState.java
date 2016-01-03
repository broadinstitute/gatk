package org.broadinstitute.hellbender.utils.hmm;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Three-state copy number hidden states enum.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public enum CopyNumberTriState {

    /**
     * Represents a loss of one or more copies in a region.
     */
    DELETION("-"),

    /**
     * Represents not a loss nor a gain of any copy in a region.
     */
    NEUTRAL("0"),

    /**
     * Represents a gain of one or more copies in a region.
     */
    DUPLICATION("+");

    private final String callString;

    CopyNumberTriState(final String callString) {
        this.callString = callString;
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
