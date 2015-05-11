package org.broadinstitute.hellbender.tools.picard.vcf.concordance;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * Class that holds summary metrics about Genotype Concordance
 *
 * @author George Grant
 */
public final class GenotypeConcordanceSummaryMetrics extends MetricBase {
    /**
     * Empty constructor - needed for unit tests
     */
    public GenotypeConcordanceSummaryMetrics() {
    }

    GenotypeConcordanceSummaryMetrics(final VariantContext.Type variantType, final GenotypeConcordanceCounts concordanceCounts,
                                      final String truthSample, final String callSample) {
        this.VARIANT_TYPE = variantType;
        this.TRUTH_SAMPLE = truthSample;
        this.CALL_SAMPLE = callSample;

        final GenotypeConcordanceScheme scheme = new GenotypeConcordanceScheme();
        concordanceCounts.validateCountsAgainstScheme(scheme);

        this.HET_SENSITIVITY = concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HET_TRUTH_STATES);
        this.HET_PPV = concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HET_CALL_STATES);
        this.HET_SPECIFICITY = concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.HET_TRUTH_STATES);

        this.HOMVAR_SENSITIVITY = concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES);
        this.HOMVAR_PPV = concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.HOM_VAR_CALL_STATES);
        this.HOMVAR_SPECIFICITY = concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.HOM_VAR_TRUTH_STATES);

        this.VAR_SENSITIVITY = concordanceCounts.getSensitivity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES);
        this.VAR_PPV = concordanceCounts.Ppv(scheme, GenotypeConcordanceCounts.VAR_CALL_STATES);
        this.VAR_SPECIFICITY = concordanceCounts.getSpecificity(scheme, GenotypeConcordanceCounts.VAR_TRUTH_STATES);
    }

    /** The type of the event (i.e. either SNP or INDEL) */
    public VariantContext.Type VARIANT_TYPE;

    /** The name of the 'truth' sample */
    public String TRUTH_SAMPLE;

    /** The name of the 'call' sample */
    public String CALL_SAMPLE;

    /** The sensitivity for all heterozygous variants (Sensitivity is TP / (TP + FN)) */
    public double HET_SENSITIVITY;

    /** The ppv (positive predictive value) for all heterozygous variants (PPV is the TP / (TP + FP)) */
    public double HET_PPV;

    /** The specificity for all heterozygous variants (Specificity is TN / (FP + TN)) */
    public double HET_SPECIFICITY;

    /** The sensitivity for all homozygous variants (Sensitivity is TP / (TP + FN)) */
    public double HOMVAR_SENSITIVITY;

    /** The ppv (positive predictive value) for all homozygous variants (PPV is the TP / (TP + FP)) */
    public double HOMVAR_PPV;

    /** The specificity for all homozygous variants (Specificity is TN / (FP + TN)) */
    public double HOMVAR_SPECIFICITY;

    /** The sensitivity for all (heterozygous and homozygous) variants (Sensitivity is TP / (TP + FN)) */
    public double VAR_SENSITIVITY;

    /** The ppv (positive predictive value) for all (heterozygous and homozygous) variants (PPV is the TP / (TP + FP)) */
    public double VAR_PPV;

    /** The specificity for all (heterozygous and homozygous) variants (Specificity is TN / (FP + TN)) */
    public double VAR_SPECIFICITY;
}
