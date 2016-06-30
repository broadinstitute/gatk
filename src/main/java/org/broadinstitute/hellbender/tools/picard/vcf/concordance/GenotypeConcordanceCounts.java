package org.broadinstitute.hellbender.tools.picard.vcf.concordance;

import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.picard.vcf.concordance.GenotypeConcordanceStates.CallState;
import org.broadinstitute.hellbender.tools.picard.vcf.concordance.GenotypeConcordanceStates.ContingencyState;
import org.broadinstitute.hellbender.tools.picard.vcf.concordance.GenotypeConcordanceStates.TruthAndCallStates;
import org.broadinstitute.hellbender.tools.picard.vcf.concordance.GenotypeConcordanceStates.TruthState;

import java.util.*;

/**
 * A class to store the counts for various truth and call state classifications relative to a reference.  With these counts and a provided
 * scheme, summary metrics can be returned.
 * @author nhomer
 */
public final class GenotypeConcordanceCounts {

    /**
     * Pre-defined sets based on if the caller wishes to return the sensitivity given the common homozygous reference, heterozygous, and homozygous variant cases.
     */
    static final TruthState[] HOM_REF_TRUTH_STATES = {TruthState.HOM_REF};
    static final TruthState[] HET_TRUTH_STATES = {TruthState.HET_REF_VAR1, TruthState.HET_VAR1_VAR2};
    static final TruthState[] HOM_VAR_TRUTH_STATES = {TruthState.HOM_VAR1};
    static final TruthState[] VAR_TRUTH_STATES = {TruthState.HET_REF_VAR1, TruthState.HET_VAR1_VAR2, TruthState.HOM_VAR1};

    /**
     * Pre-defined sets based on if the caller wishes to return the PPV given the common homozygous reference, heterozygous, and homozygous variant cases.
     */
    static final CallState[] HOM_REF_CALL_STATES = {CallState.HOM_REF};
    static final CallState[] HET_CALL_STATES = {CallState.HET_REF_VAR1, CallState.HET_REF_VAR2, CallState.HET_REF_VAR3,
            CallState.HET_VAR1_VAR2, CallState.HET_VAR1_VAR3, CallState.HET_VAR3_VAR4};
    static final CallState[] HOM_VAR_CALL_STATES = {CallState.HOM_VAR1, CallState.HOM_VAR2, CallState.HOM_VAR3};
    static final CallState[] VAR_CALL_STATES = {CallState.HET_REF_VAR1, CallState.HET_REF_VAR2, CallState.HET_REF_VAR3,
            CallState.HET_VAR1_VAR2, CallState.HET_VAR1_VAR3, CallState.HET_VAR3_VAR4,
            CallState.HOM_VAR1, CallState.HOM_VAR2, CallState.HOM_VAR3};

    /** The underlying counts table */
    private final Histogram<TruthAndCallStates> counter = new Histogram<>();

    /**
     * Increments a count for the truth/call state tuple.
     * @param truthAndCallStates
     */
    public void increment(final TruthAndCallStates truthAndCallStates) {
        this.counter.increment(truthAndCallStates);
    }

    /**
     * Validates that there are no counts for NA states in the underlying scheme
     */
    public void validateCountsAgainstScheme(final GenotypeConcordanceScheme scheme) {
        final Set<ContingencyState> naContingencyStates = getContingencyStateSet(GenotypeConcordanceScheme.NA);
        for (final TruthState truthState : TruthState.values()) {
            for (final CallState callState : CallState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                if (0 < getCount(truthAndCallStates)) {
                    final Set<ContingencyState> contingencyStates = getContingencyStateSet(scheme.getConcordanceStateArray(truthAndCallStates));
                    if (contingencyStates.containsAll(naContingencyStates)) {
                        throw new GATKException(String.format("Found counts for an illegal set of states: [%s, %s]", truthState.name(), callState.name()));
                    }
                }
            }
        }
    }

    private Set<ContingencyState> getContingencyStateSet(final ContingencyState[] contingencyStateArray) {
        final Set<ContingencyState> contingencyStateSet = new HashSet<>();
        Collections.addAll(contingencyStateSet, contingencyStateArray);
        return contingencyStateSet;
    }

    /**
     * Returns the sensitivity defined by the scheme across the subset of truth states.
     */
    public double getSensitivity(final GenotypeConcordanceScheme scheme, final TruthState[] truthStateArray) {
        /**
         * Sensitivity is the TP / P = TP / (TP + FN)
         */
        double numerator = 0.0;
        double denominator = 0.0;

        scheme.validateScheme();

        for (final TruthState truthState : truthStateArray) {
            for (final CallState callState : CallState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                final int count = getCount(truthAndCallStates);
                for (final ContingencyState contingencyState : scheme.getConcordanceStateArray(truthAndCallStates)) {
                    if (ContingencyState.TP == contingencyState) {
                        numerator += count;
                        denominator += count;
                    } else if (ContingencyState.FN == contingencyState) {
                        denominator += count;
                    }
                }
            }
        }

        return (numerator / denominator);
    }

    /**
     * Returns the PPV defined by the scheme across the subset of call states.
     */
    public double Ppv(final GenotypeConcordanceScheme scheme, final CallState[] callStateList) {
        /**
         * PPV is the TP / (TP + FP)
         */
        double numerator = 0.0;
        double denominator = 0.0;

        scheme.validateScheme();

        for (final CallState callState : callStateList) {
            for (final TruthState truthState : TruthState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                final int count = getCount(truthAndCallStates);
                for (final ContingencyState contingencyState : scheme.getConcordanceStateArray(truthAndCallStates)) {
                    if (ContingencyState.TP == contingencyState) {
                        numerator += count;
                        denominator += count;
                    } else if (ContingencyState.FP == contingencyState) {
                        denominator += count;
                    }
                }
            }
        }

        return (numerator / denominator);
    }

    /**
     * Returns the specificity defined by the scheme across the subset of truth states.
     */
    public double getSpecificity(final GenotypeConcordanceScheme scheme, final TruthState[] truthStateArray) {
        /**
         * Specificity is the TN / N = TN / (FP + TN)
         */
        double numerator = 0.0;
        double denominator = 0.0;

        scheme.validateScheme();

        for (final TruthState truthState : truthStateArray) {
            for (final CallState callState : CallState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                final int count = getCount(truthAndCallStates);
                for (final ContingencyState contingencyState : scheme.getConcordanceStateArray(truthAndCallStates)) {
                    if (ContingencyState.TN == contingencyState) {
                        numerator += count;
                        denominator += count;
                    } else if (ContingencyState.FP == contingencyState) {
                        denominator += count;
                    }
                }
            }
        }
        return (numerator / denominator);
    }

    /**
     * Returns the count defined by the truth state set and call state set.
     */
    public int getCount(final TruthState truthState, final CallState callState) {
        return getCount(new TruthAndCallStates(truthState, callState));
    }

    /**
     * Returns the count defined by the truth state set and call state set.
     */

    @SuppressWarnings("unchecked")
    public int getCount(final TruthAndCallStates truthAndCallStates) {
        final Histogram.Bin<TruthAndCallStates> bin = this.counter.get(truthAndCallStates);
        return (bin == null ? 0 : (int) bin.getValue());
    }



    /**
     * Returns the sum of all pairs of tuples defined by the truth state set and call state set.
     */
    public int getSum(final Set<TruthState> truthStateSet, final Set<CallState> callStateSet) {
        int count = 0;
        for (final TruthState truthState : truthStateSet) {
            for (final CallState callState : callStateSet) {
                count += getCount(truthState, callState);
            }
        }
        return count;
    }

    /**
     * Returns the total number of times each contingency state is encountered, summed across all truth/call state pairs.
     */
    public Map<ContingencyState, Integer> getContingencyStateCounts(final GenotypeConcordanceScheme scheme) {
        scheme.validateScheme();

        final Map<ContingencyState, Integer> counts = new HashMap<>();
        for (final ContingencyState contingencyState : ContingencyState.values()) {
            counts.put(contingencyState, 0);
        }

        for (final TruthState truthState : TruthState.values()) {
            for (final CallState callState : CallState.values()) {
                final TruthAndCallStates truthAndCallStates = new TruthAndCallStates(truthState, callState);
                final ContingencyState[] contingencyStateArray = scheme.getConcordanceStateArray(truthAndCallStates);
                for (final ContingencyState contingencyState : contingencyStateArray) {
                    final int newCount = counts.get(contingencyState) + getCount(truthAndCallStates);
                    counts.put(contingencyState, newCount);
                }
            }
        }

        return counts;
    }
}
