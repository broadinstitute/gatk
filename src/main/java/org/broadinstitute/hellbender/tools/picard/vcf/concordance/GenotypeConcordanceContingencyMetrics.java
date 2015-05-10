package org.broadinstitute.hellbender.tools.picard.vcf.concordance;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.picard.vcf.concordance.GenotypeConcordanceStates.*;

import java.util.*;

/**
 * Class that holds metrics about the Genotype Concordance contingency tables.
 */
public final class GenotypeConcordanceContingencyMetrics extends MetricBase {
    /**
     * Empty constructor - needed for unit tests
     */
    public GenotypeConcordanceContingencyMetrics() {
    }

    GenotypeConcordanceContingencyMetrics(final VariantContext.Type variantType, final GenotypeConcordanceCounts concordanceCounts,
                                          final String truthSample, final String callSample) {
        this.VARIANT_TYPE = variantType;
        this.TRUTH_SAMPLE = truthSample;
        this.CALL_SAMPLE = callSample;

        final GenotypeConcordanceScheme scheme = new GenotypeConcordanceScheme();
        concordanceCounts.validateCountsAgainstScheme(scheme);

        Map<ContingencyState, Integer> counts = concordanceCounts.getContingencyStateCounts(scheme);
        this.TP_COUNT = counts.get(ContingencyState.TP);
        this.TN_COUNT = counts.get(ContingencyState.TN);
        this.FP_COUNT = counts.get(ContingencyState.FP);
        this.FN_COUNT = counts.get(ContingencyState.FN);
        this.EMPTY_COUNT = counts.get(ContingencyState.EMPTY);
    }

    /** The type of the event (i.e. either SNP or INDEL) */
    public VariantContext.Type VARIANT_TYPE;

    /** The name of the 'truth' sample */
    public String TRUTH_SAMPLE;

    /** The name of the 'call' sample */
    public String CALL_SAMPLE;

    /** The TP (true positive) count across all variants */
    public int TP_COUNT;

    /** The TN (true negative) count across all variants */
    public int TN_COUNT;

    /** The FP (false positive) count across all variants */
    public int FP_COUNT;

    /** The FN (false negative) count across all variants */
    public int FN_COUNT;

    /** The empty (no contingency info) count across all variants */
    public int EMPTY_COUNT;
}
