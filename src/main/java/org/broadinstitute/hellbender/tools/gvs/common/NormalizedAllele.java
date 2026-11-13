package org.broadinstitute.hellbender.tools.gvs.common;

/**
 * Result of VRS allele normalization.
 * 
 * Represents a normalized allele with its location coordinates and state.
 * This is an intermediate representation before computing the final VRS ID.
 */
public class NormalizedAllele {
    public final long start;  // 0-based interbase
    public final long end;    // 0-based interbase
    public final SequenceState state;

    public NormalizedAllele(long start, long end, SequenceState state) {
        this.start = start;
        this.end = end;
        this.state = state;
    }

    /**
     * Represents the state of an allele, which can be either a literal sequence
     * or a reference length expression.
     */
    public static class SequenceState {
        public final StateType type;
        public final String sequence;        // For LITERAL_SEQUENCE_EXPRESSION
        public final Integer length;         // For REFERENCE_LENGTH_EXPRESSION
        public final Integer repeatSubunitLength;  // For REFERENCE_LENGTH_EXPRESSION

        private SequenceState(StateType type, String sequence, Integer length, Integer repeatSubunitLength) {
            this.type = type;
            this.sequence = sequence;
            this.length = length;
            this.repeatSubunitLength = repeatSubunitLength;
        }

        /**
         * Create a literal sequence state (for SNPs and simple indels).
         */
        public static SequenceState literalSequence(String sequence) {
            return new SequenceState(StateType.LITERAL_SEQUENCE_EXPRESSION, sequence, null, null);
        }

        /**
         * Create a reference length state (for tandem repeat expansions).
         */
        public static SequenceState referenceLength(int length, int repeatSubunitLength) {
            return new SequenceState(StateType.REFERENCE_LENGTH_EXPRESSION, null, length, repeatSubunitLength);
        }

        public enum StateType {
            LITERAL_SEQUENCE_EXPRESSION,
            REFERENCE_LENGTH_EXPRESSION
        }
    }
}
