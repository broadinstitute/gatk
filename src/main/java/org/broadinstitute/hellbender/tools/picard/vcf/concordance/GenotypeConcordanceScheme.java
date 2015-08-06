package org.broadinstitute.hellbender.tools.picard.vcf.concordance;

import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.picard.vcf.concordance.GenotypeConcordanceStates.*;

import java.util.*;

/**
 * This defines for each valid TruthState and CallState tuple, the set of contingency table entries that to which the tuple should contribute.
 * @author nhomer
 */
public final class GenotypeConcordanceScheme {

    /** The underlying scheme */
    protected final Map<TruthAndCallStates, ContingencyState[]> scheme = new HashMap<>();

    /** These are convenience variables for defining a scheme.  NA means that such a tuple should never be observed. */
    public static final ContingencyState[]    NA       = {ContingencyState.NA};
    protected static final ContingencyState[] EMPTY    = {ContingencyState.EMPTY};
    protected static final ContingencyState[] TP_ONLY  = {ContingencyState.TP};
    protected static final ContingencyState[] FP_ONLY  = {ContingencyState.FP};
    protected static final ContingencyState[] TN_ONLY  = {ContingencyState.TN};
    protected static final ContingencyState[] FN_ONLY  = {ContingencyState.FN};
    protected static final ContingencyState[] TP_FN    = {ContingencyState.TP, ContingencyState.FN};
    protected static final ContingencyState[] TP_FP    = {ContingencyState.TP, ContingencyState.FP};
    protected static final ContingencyState[] TP_TN    = {ContingencyState.TP, ContingencyState.TN};
    protected static final ContingencyState[] FP_FN    = {ContingencyState.FP, ContingencyState.FN};
    protected static final ContingencyState[] FP_TN    = {ContingencyState.FP, ContingencyState.TN};
    protected static final ContingencyState[] FP_TN_FN = {ContingencyState.FP, ContingencyState.TN, ContingencyState.FN};
    protected static final ContingencyState[] TP_FP_FN = {ContingencyState.TP, ContingencyState.FP, ContingencyState.FN};
    protected static final ContingencyState[] TN_FN    = {ContingencyState.TN, ContingencyState.FN};


    /** Has this scheme been previously validated */
    private boolean isValidated = false;

    /**
     * Adds a row to the scheme
     * @param callState the call state (row)
     * @param concordanceStateArrays the concordance state arrays for each truth value, in order
     */
    protected void addRow(final CallState callState, final ContingencyState[]... concordanceStateArrays) {
        if (concordanceStateArrays.length != GenotypeConcordanceStates.TruthState.values().length) {
            throw new GATKException("Length mismatch between concordanceStateArrays and TruthState.values()");
        }
        for (int i = 0; i < concordanceStateArrays.length; i++) {
            scheme.put(new TruthAndCallStates(TruthState.values()[i], callState), concordanceStateArrays[i]);
        }
    }

    /**
     * The scheme is defined in the constructor.
     *
     * The default scheme is derived from the GA4GH Benchmarking Work Group's proposed evaluation scheme.
     *
     * In general, we are comparing two sets of alleles.  Therefore, we can have zero or more contingency table values represented in one comparison.  For example, if the truthset is
     * a heterozygous call with both alleles non-reference (HET_VAR1_VAR2), and the callset is a heterozygous call with both alleles non-reference with one of the alternate alleles
     * matching an alternate allele in the callset, we would have a true positive, false positive, and false negative.  The true positive is from the matching alternate alleles, the
     * false positive is the alternate allele found in the callset but not found in the truthset, and the false negative is the alternate in the truthset not found in the callset.
     *
     * We also include a true negative in cases where the reference allele is found in both the truthset and callset.
     *
     * We have no HET_VAR2_VAR3 case, as VAR2/VAR3 are simply symbolic, and so we can change HET_VAR2_VAR3 into the HET_VAR3_VAR4 case.
     *
     * In this (the default) scheme
     *
     * Finally, we have NA cases, which represent tuples that our code can and should not reach.
     */
    public GenotypeConcordanceScheme() {

        /**          ROW STATE            MISSING       HOM_REF       HET_REF_VAR1       HET_VAR1_VAR2        HOM_VAR1        NO_CALL        LOW_GQ        LOW_DP        VC_FILTERED   GT_FILTERED   IS_MIXED    **/
        addRow(CallState.MISSING,         NA,           TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HOM_REF,         TN_ONLY,      TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_REF_VAR1,    FP_TN,        FP_TN,        TP_TN,             TP_FN,               TP_FN,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_REF_VAR2,    NA,           NA,           FP_TN_FN,          NA,                  FP_FN,          NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.HET_REF_VAR3,    NA,           NA,           NA,                FP_FN,               NA,             NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.HET_VAR1_VAR2,   FP_ONLY,      FP_ONLY,      TP_FP,             TP_ONLY,             TP_FP_FN,       EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HET_VAR1_VAR3,   NA,           NA,           NA,                TP_FP_FN,            NA,             NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.HET_VAR3_VAR4,   FP_ONLY,      FP_ONLY,      FP_FN,             FP_FN,               FP_FN,          NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.HOM_VAR1,        FP_ONLY,      FP_ONLY,      TP_FP,             TP_FN,               TP_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.HOM_VAR2,        NA,           NA,           FP_FN,             TP_FN,               FP_FN,          NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.HOM_VAR3,        NA,           NA,           NA,                FP_FN,               NA,             NA,            NA,           NA,           NA,           NA,           NA);
        addRow(CallState.NO_CALL,         EMPTY,        EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.VC_FILTERED,     EMPTY,        TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.GT_FILTERED,     EMPTY,        TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.LOW_GQ,          EMPTY,        TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.LOW_DP,          EMPTY,        TN_ONLY,      TN_FN,             FN_ONLY,             FN_ONLY,        EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);
        addRow(CallState.IS_MIXED,        EMPTY,        EMPTY,        EMPTY,             EMPTY,               EMPTY,          EMPTY,         EMPTY,        EMPTY,        EMPTY,        EMPTY,        EMPTY);

        validateScheme();
    }

    /**
     * Get the concordance state array associate with the given truth state and call state tuple.
     */
    public ContingencyState[] getConcordanceStateArray(final TruthState truthState, final CallState callState) {
        return this.getConcordanceStateArray(new TruthAndCallStates(truthState, callState));
    }

    /**
     * Get the concordance state array associate with the given truth state and call state tuple.
     */
    public ContingencyState[] getConcordanceStateArray(final TruthAndCallStates truthAndCallStates) {
        return this.scheme.get(truthAndCallStates);
    }

    /**
     * Get the contingency state array as a parse-able string
     */
    public String getContingencyStateString(final TruthState truthState, final CallState callState) {
        final ContingencyState[] contingencyStateArray = getConcordanceStateArray(truthState, callState);
        return (contingencyStateArray.length == 0) ? "EMPTY" : StringUtil.join(",", contingencyStateArray);
    }

    /**
     * Check that all cells in the scheme exist.
     * @throws GATKException if a missing tuple was found.
     */
    public void validateScheme() throws GATKException {
        if (!isValidated) {
            for (final TruthState truthState : TruthState.values()) {
                for (final CallState callState : CallState.values()) {
                    if (!scheme.containsKey(new TruthAndCallStates(truthState, callState))) {
                        throw new GATKException(String.format("Missing scheme tuple: [%s, %s]", truthState.name(), callState.name()));
                    }
                }
            }
        }

        isValidated = true;
    }
}
