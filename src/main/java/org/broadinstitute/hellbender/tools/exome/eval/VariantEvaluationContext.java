package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriStateAllele;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.GenotypeCopyNumberTriStateSegments;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

/**
 * Enhanced variant context for copy-number-tri state variants with evaluation information.
 */
@SuppressWarnings("serial")
public class VariantEvaluationContext extends VariantContext {

    public static final String TRUTH_ALLELE_FREQUENCY_KEY = "TruthAF";
    public static final String CALLS_ALLELE_FREQUENCY_KEY = VCFConstants.ALLELE_FREQUENCY_KEY;
    public static final String TRUTH_GENOTYPE_KEY = "TGT";
    public static final String CALLS_ALLELE_NUMBER_KEY = VCFConstants.ALLELE_NUMBER_KEY;
    public static final String TRUTH_ALLELE_NUMBER_KEY = "TruthAN";
    public static final String EVALUATION_CLASS_KEY = "EV";
    public static final String CALLED_SEGMENTS_COUNT_KEY = "SC";
    public static final String CALLED_ALLELE_COUNTS_KEY = "CC";
    public static final String CALLED_TARGET_COUNT_KEY = "CL";
    public static final String CALL_QUALITY_KEY = "CQ";
    public static final String TRUTH_TARGET_COUNT_KEY = GenotypeCopyNumberTriStateSegments.NUMBER_OF_TARGETS_KEY;
    public static final String TRUTH_COPY_FRACTION_KEY = "TF";
    public static final String TRUTH_QUALITY_KEY = "TQ";

    private final int truthAN;
    private final double[] truthAF;
    private final double[] callsAF;
    private final int targetCount;
    private VariantContext truth;
    private List<VariantContext> calls;

    protected VariantEvaluationContext(final VariantContext other) {
        super(other);
        truthAN = getAttributeAsInt(TRUTH_ALLELE_NUMBER_KEY, 0);
        truthAF = getAlleleDoubleArrayFromAlternativeAlleleArrayAttribute(TRUTH_ALLELE_FREQUENCY_KEY, Double.NaN);
        callsAF = getAlleleDoubleArrayFromAlternativeAlleleArrayAttribute(CALLS_ALLELE_FREQUENCY_KEY, Double.NaN);
        targetCount = getAttributeAsInt(TRUTH_TARGET_COUNT_KEY, 0);
    }

    public int getTargetCount() {
        return targetCount;
    }

    public int getTruthAlleleNumber() {
        return truthAN;
    }

    public double getTruthAlleleFrequency(final CopyNumberTriStateAllele allele) {
        Utils.nonNull(allele);
        final int index = alleles.indexOf(allele);
        return index >= 0 ? truthAF[index] : 0;
    }

    public double getCallsAlleleFrequency(final CopyNumberTriStateAllele allele) {
        Utils.nonNull(allele);
        final int index = alleles.indexOf(allele);
        return index >= 0 ? callsAF[index] : 0;
    }

    private double[] getAlleleDoubleArrayFromAlternativeAlleleArrayAttribute(final String key, final double missingValue) {
        final double[] alternativeAllelesFrequencies = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(this, key,
                () -> { return new double[alleles.size() - 1]; }, missingValue);
        final double[] result = new double[alleles.size()];
        if (alternativeAllelesFrequencies.length != alleles.size() - 1) {
            throw new IllegalStateException(String.format("We expect the %s Info annotation to contain an array of %d elements (alt. allele count)", key, alleles.size() - 1));
        }
        System.arraycopy(alternativeAllelesFrequencies, 0, result, 1, alternativeAllelesFrequencies.length);
        final double nonRefSum = MathUtils.sum(result);
        Utils.validateArg(nonRefSum <= 1.0, () -> String.format("The sum of element on annotation %s cannot greater than 1.0: %g", key, nonRefSum));
        result[0] = 1.0 - nonRefSum;
        return result;
    }

    protected void setEvidence(final VariantContext truth, final List<VariantContext> calls) {
        this.truth = truth;
        this.calls = calls;
    }

    public VariantContext getTruthVariantContext() {
        return truth;
    }

    public List<VariantContext> getCallsVariantContexts() {
        return calls;
    }
}
