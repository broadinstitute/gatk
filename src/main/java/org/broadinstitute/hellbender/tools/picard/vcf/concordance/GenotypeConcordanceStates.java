package org.broadinstitute.hellbender.tools.picard.vcf.concordance;

/**
 * A class to store the various classifications for:
 * 1. a truth genotype versus a reference
 * 2. a call genotype versus a truth, relative to a reference
 *
 * An example use of this class is to have one instance per following use case:
 * - SNPs
 * - indels
 * - filtered variant (truth or call)
 * - filtered genotype (truth or call)
 * - low GQ (call)
 * - low DP (call)
 * - No call (truth or call)
 * - No variant (truth or call) *
 *
 * @author nhomer
 */
public final class GenotypeConcordanceStates {
    /**
     * These states represent the relationship between a truth genotype and the reference sequence.
     */
    public enum TruthState {
        MISSING,
        HOM_REF, // ref/ref
        HET_REF_VAR1, // ref/var1 (var1!=ref)
        HET_VAR1_VAR2, // var1/var2 (var1!=var2, var1!=ref, var2!=ref)
        HOM_VAR1, // var1/var1 (var1!=ref)
        NO_CALL,
        LOW_GQ,
        LOW_DP,
        VC_FILTERED,
        GT_FILTERED,
        IS_MIXED;

        public static TruthState getHom(final int alleleIdx) {
            if (alleleIdx == 0) return HOM_REF;
            if (alleleIdx == 1) return HOM_VAR1;
            assert false;
            return null;
        }

        public static TruthState getVar(final int allele0idx, final int allele1idx) {
            if (allele0idx == 0 && allele1idx == 1) return HET_REF_VAR1;
            if (allele0idx == 1 && allele1idx == 0) return HET_REF_VAR1;

            if (allele0idx == 1 && allele1idx == 2) return HET_VAR1_VAR2;
            if (allele0idx == 2 && allele1idx == 1) return HET_VAR1_VAR2;

            assert false;
            return null;
        }
    }

    /**
     * These states represent the relationship between the call genotype and the truth genotype relative to
     * a reference sequence.
     */
    enum CallState {
        MISSING,
        HOM_REF, // ref/ref, valid for all TruthStates
        HET_REF_VAR1, // ref/var1, valid for all TruthStates
        HET_REF_VAR2, // ref/var2, valid only for TruthStates: HET_REF_VAR1, HET_VAR1_VAR2, HOM_VAR1
        HET_REF_VAR3, // ref/var3, valid only for TruthStates: HET_VAR1_VAR2
        HET_VAR1_VAR2, // var1/var2, valid for all TruthStates
        HET_VAR1_VAR3, // var1/var3, valid only for TruthStates: HET_VAR1_VAR2. also encapsulates HET_VAR2_VAR3 (see special case below)
        HET_VAR3_VAR4, // var3/var4, valid only for TruthStates: HET_REF_VAR1, HET_VAR1_VAR2, HOM_VAR1
        HOM_VAR1, // var1/var1, valid for all TruthStates
        HOM_VAR2, // var2/var2, valid only for TruthStates: HET_REF_VAR1, HET_VAR1_VAR2, HOM_VAR1
        HOM_VAR3, // var3/var3, valid only for TruthStates: HET_VAR1_VAR2
        NO_CALL,
        LOW_GQ,
        LOW_DP,
        VC_FILTERED,
        GT_FILTERED,
        IS_MIXED;

        public static CallState getHom(final int alleleIdx) {
            if (alleleIdx == 0) return HOM_REF;
            if (alleleIdx == 1) return HOM_VAR1;
            if (alleleIdx == 2) return HOM_VAR2;
            if (alleleIdx == 3) return HOM_VAR3;

            assert false;
            return null;
        }

        public static CallState getHet(int allele0idx, int allele1idx) {

            if(allele0idx > allele1idx){
                final int temp = allele0idx;
                allele0idx=allele1idx;
                allele1idx=temp;
            }
            if(allele0idx == 0) { //REF CASE
                if (allele1idx == 1) return HET_REF_VAR1;
                if (allele1idx == 2) return HET_REF_VAR2;
                if (allele1idx == 3) return HET_REF_VAR3;
                assert false;
                return null;
            }

            //HET CASES
            if(allele0idx == 1) {
                if (allele1idx == 2) return HET_VAR1_VAR2;
                if (allele1idx == 3) return HET_VAR1_VAR3;
                assert false;
                return null;
            }

            if(allele0idx == 2 && allele1idx == 3) return HET_VAR3_VAR4; //special case not a mistake.
            if(allele0idx == 3 && allele1idx == 4) return HET_VAR3_VAR4;

            assert false;
            return null;
        }
    }

    /**
     * A specific state for a 2x2 contingency table.
     * NA denotes an invalid state that should not be reachable by the code.
     * EMPTY denotes that no conclusion could be drawn from the data.
     */
    enum ContingencyState {
        TP,
        FP,
        TN,
        FN,
        NA,
        EMPTY
    }

    /**
     * A minute class to store the truth and call state respectively.
     */
    static class TruthAndCallStates implements Comparable<TruthAndCallStates> {
        public final TruthState truthState;
        public final CallState callState;

        public TruthAndCallStates(final TruthState truthState, final CallState callState) {
            this.truthState = truthState;
            this.callState = callState;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            return compareTo((TruthAndCallStates) o) == 0;
        }

        @Override
        public int hashCode() {
            int result = truthState.hashCode();
            result = 31 * result + callState.hashCode();
            return result;
        }

        @Override
        public int compareTo(final TruthAndCallStates that) {
            int result = this.truthState.compareTo(that.truthState);
            if (result == 0) result = this.callState.compareTo(that.callState);
            return result;
        }
    }
}
