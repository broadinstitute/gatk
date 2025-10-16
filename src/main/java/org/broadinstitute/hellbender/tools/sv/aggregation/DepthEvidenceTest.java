package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.DepthEvidenceGenotyper;
import org.broadinstitute.hellbender.tools.sv.DepthMatrix;
import org.broadinstitute.hellbender.tools.sv.PermutationTTest;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.hipparchus.stat.descriptive.moment.StandardDeviation;
import picard.util.MathUtil;

import java.util.*;
import java.util.stream.Collectors;

public class DepthEvidenceTest {

    protected static final Median MEDIAN = new Median();
    protected static final NormalDistribution STD_NORMAL = new NormalDistribution();
    protected static final double POWER_SIGNIFICANCE_LEVEL = 0.05;

    protected final double powerThreshold;

    protected record SampleSetResult(Set<String> treatSampleSet, Set<String> controlSampleSet) {
    }

    protected record DepthStats(double pValue, double secondMaxP) {
    }

    public record DepthTestResult(double pValue, double secondMaxP, double medianSeparation) {
    }

    public DepthEvidenceTest(final double powerThreshold) {
        this.powerThreshold = powerThreshold;
    }

    public DepthTestResult test(final SVCallRecord record, final DepthMatrix depthMatrix) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType svType = record.getType();
        Utils.validateArg(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL.equals(svType) || GATKSVVCFConstants.StructuralVariantAnnotationType.DUP.equals(svType), "Only DEL/DUP types are supported");
        final SampleSetResult sampleSets = getSampleSets(record);
        if (sampleSets == null) return null;
        final double[] controlMedians = getSampleMedians(depthMatrix, sampleSets.controlSampleSet);
        final double[] treatMedians = getSampleMedians(depthMatrix, sampleSets.treatSampleSet);
        final double power = power(controlMedians, treatMedians, POWER_SIGNIFICANCE_LEVEL);
        final DepthStats result;
        if (!Double.isNaN(power) && power > powerThreshold) {
            result = poweredTest(depthMatrix, sampleSets.controlSampleSet, sampleSets.treatSampleSet, svType, controlMedians, treatMedians);
        } else {
            result = underpoweredTest(depthMatrix, sampleSets.controlSampleSet, sampleSets.treatSampleSet, svType, treatMedians, controlMedians);
        }
        final double medianSeparation = getMedianSeparation(svType, controlMedians, treatMedians);
        return new DepthTestResult(result.pValue, result.secondMaxP, medianSeparation);
    }

    public SVCallRecord applyToRecord(final SVCallRecord record,
                                      final DepthTestResult result,
                                      final double maxQual,
                                      final SAMSequenceDictionary dictionary) {
        Utils.nonNull(record);
        Utils.nonNull(result);
        final Map<String, Object> refinedAttr = new HashMap<>(record.getAttributes());
        refinedAttr.put(GATKSVVCFConstants.READ_DEPTH_QUALITY_ATTRIBUTE, Math.min(-10. * Math.log10(result.pValue()), maxQual));
        refinedAttr.put(GATKSVVCFConstants.READ_DEPTH_SECOND_MAX_QUALITY_ATTRIBUTE, Math.min(-10. * Math.log10(result.secondMaxP()), maxQual));
        refinedAttr.put(GATKSVVCFConstants.READ_DEPTH_MEDIAN_SEPARATION_ATTRIBUTE, result.medianSeparation());

        // Create new record
        return new SVCallRecord(record.getId(), record.getContigA(), record.getPositionA(),
                record.getStrandA(), record.getContigB(), record.getPositionB(), record.getStrandB(),
                record.getType(), record.getComplexSubtype(), record.getComplexEventIntervals(), record.getLength(),
                record.getEvidence(), record.getAlgorithms(), record.getAlleles(),
                record.getGenotypes(), refinedAttr, record.getFilters(), record.getLog10PError(), dictionary);
    }

    protected static SampleSetResult getSampleSets(final SVCallRecord record) {
        final List<Genotype> carrierGenotypes = record.getCarrierGenotypeList();
        if (carrierGenotypes.isEmpty()) {
            // Do not process uncalled variants
            return null;
        }
        final int maxExpectedCopyNumber = getMaxExpectedCopyNumber(carrierGenotypes);
        if (maxExpectedCopyNumber > 2) {
            throw new GATKException("Found a sample with expected copy number (ECN) of " + maxExpectedCopyNumber + " but this tool only supports diploid samples");
        }
        final Set<String> rawCarrierSampleSet = record.getCarrierSampleSet();
        final Set<String> rawRefSampleSet = Sets.difference(record.getAllSamples(), rawCarrierSampleSet); // TODO: assumes no-calls are hom-ref
        final Set<String> eligibleGenotypes = getSamplesWithExpectedCopyNumber(maxExpectedCopyNumber, record.getGenotypes()); // TODO: male-only variants only tested against other males but females could be incorporated here
        final Set<String> treatSampleSet = Sets.intersection(rawCarrierSampleSet, eligibleGenotypes);
        final Set<String> controlSampleSet = Sets.intersection(rawRefSampleSet, eligibleGenotypes);
        return new SampleSetResult(treatSampleSet, controlSampleSet);
    }

    private static int getMaxExpectedCopyNumber(final Collection<Genotype> genotypes) {
        Utils.nonNull(genotypes);
        if (genotypes.isEmpty()) {
            return 0;
        }
        return genotypes.stream()
                .mapToInt(g -> VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1))
                .max().getAsInt();
    }

    private static Set<String> getSamplesWithExpectedCopyNumber(final int state, final List<Genotype> genotypes) {
        Utils.nonNull(genotypes);
        if (genotypes.isEmpty()) {
            return Collections.emptySet();
        }
        return genotypes.stream()
                .filter(g -> VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, -1) == state)
                .map(Genotype::getSampleName)
                .collect(Collectors.toUnmodifiableSet());
    }

    public static double[] getSampleMedians(final DepthMatrix depthMatrix, final Collection<String> samples) {
        final double[] medians = new double[samples.size()];
        int i = 0;
        for (final String sample : samples) {
            medians[i++] = MEDIAN.evaluate(depthMatrix.getSample(sample));
        }
        return medians;
    }

    protected static double power(double[] control, double[] treat, final double sigLevel) {
        Utils.nonNull(control);
        Utils.nonNull(treat);
        Utils.validateArg(sigLevel > 0 && sigLevel < 1, "significance level not on (0, 1)");

        // Set treatment group to be the smaller
        if (treat.length > control.length) {
            final double[] tmp = control;
            control = treat;
            treat = tmp;
        }

        // Calculate standard deviation of control group
        final double controlSd = new StandardDeviation().evaluate(control) * (control.length - 1.0) / control.length;

        // Avoid division by zero
        if (controlSd == 0.0 || Double.isNaN(controlSd)) {
            return Double.NaN;
        }

        // Calculate effect size (Cohen's d)
        final double effectSize = 0.5 / controlSd;

        // Calculate power using two-sample t-test with unequal sample sizes
        return powerT2nTest(control.length, treat.length, sigLevel, effectSize);
    }

    /**
     * Calculate statistical power for two-sample t-test with unequal sample sizes
     * This is equivalent to R's pwr.t2n.test function, assuming "greater" alternative hypothesis
     *
     * @param n1 Sample size of group 1
     * @param n2 Sample size of group 2
     * @param sigLevel Significance level (alpha)
     * @param effectSize Effect size (Cohen's d)
     * @return Statistical power
     */
    private static double powerT2nTest(int n1, int n2, double sigLevel, double effectSize) {
        Utils.validateArg(sigLevel > 0 && sigLevel < 1, "significance level not on (0, 1)");
        Utils.validateArg(effectSize > 0, "effect size must be positive");

        // Not enough samples
        if (n1 < 2 || n2 < 2) {
            return Double.NaN;
        }

        // Calculate degrees of freedom
        final double df = n1 + n2 - 2;

        // Calculate pooled standard error
        final double se = Math.sqrt((1.0 / n1) + (1.0 / n2));

        // Calculate non-centrality parameter
        final double ncp = effectSize / se;

        // Calculate critical value based on alternative hypothesis and significance level
        // Assume "greater than" alternative hypothesis
        final TDistribution tDistribution = new TDistribution(df);
        final double criticalValue = tDistribution.inverseCumulativeProbability(1 - sigLevel);

        // Calculate power using non-central t-distribution
        final double power = 1 - approximateNonCentralTCdf(criticalValue, df, ncp);

        return Math.max(0.0, Math.min(1.0, power)); // Ensure power is between 0 and 1
    }

    /**
     * Non-central t-distribution cumulative distribution function
     * Approximation for moderate non-centrality parameters
     *
     * @param t t-value
     * @param df Degrees of freedom
     * @param ncp Non-centrality parameter
     * @return Cumulative probability
     */
    private static double approximateNonCentralTCdf(double t, double df, double ncp) {
        // Simple approximation using shifted central t-distribution
        final double shiftedT = t - ncp * Math.sqrt(df / (df + t * t));
        return new TDistribution(df).cumulativeProbability(shiftedT);
    }

    protected static DepthStats underpoweredTest(final DepthMatrix depthMatrix,
                                                 final Set<String> controlSampleSet,
                                                 final Set<String> treatSampleSet,
                                                 final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                                                 final double[] treatMedians,
                                                 final double[] controlMedians) {
        // Underpowered case - use multiple single-sample t-tests
        final int alternativeSign;
        if (svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
            alternativeSign = 1;
        } else if (svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            alternativeSign = -1;
        } else {
            throw new IllegalStateException("Unsupported variant type: " + svtype);
        }
        final double[] pValues = new double[treatMedians.length];
        final double[] secondMaxPValues = new double[treatMedians.length];
        final double controlMean = MathUtil.mean(controlMedians);
        final double controlStd = MathUtil.stddev(controlMedians, controlMean) * Math.sqrt(controlMedians.length / (double) (controlMedians.length - 1));
        final String[] treatSampleList = treatSampleSet.toArray(new String[0]);
        final int numBins = depthMatrix.getNumBins();
        for (int i = 0; i < treatMedians.length; i++) {
            pValues[i] = Math.log(Math.max(STD_NORMAL.cumulativeProbability(alternativeSign * (treatMedians[i] - controlMean) / controlStd), Double.MIN_VALUE));
            final double[] slicePvals = new double[numBins];
            for (int j = 0; j < numBins; j++) {
                final double[] controlSlice = depthMatrix.slice(j, controlSampleSet);
                final double[] treatSlice = depthMatrix.slice(j, Collections.singleton(treatSampleList[i]));
                final double controlMeanSlice = MathUtil.mean(controlSlice);
                final double controlStdSlice = MathUtil.stddev(controlSlice, controlMeanSlice) * Math.sqrt(controlSlice.length / (double) (controlSlice.length - 1));
                slicePvals[j] = Math.log(Math.max(STD_NORMAL.cumulativeProbability(alternativeSign * (treatSlice[0] - controlMeanSlice) / controlStdSlice), Double.MIN_VALUE));
            }
            Arrays.sort(slicePvals);
            if (slicePvals.length > 0) {
                secondMaxPValues[i] = slicePvals[Math.max(slicePvals.length - 2, 0)];
            } else {
                secondMaxPValues[i] = Double.NaN;
            }
        }
        return new DepthStats(combinePValues(pValues), combinePValues(secondMaxPValues));
    }

    protected static DepthStats poweredTest(final DepthMatrix depthMatrix,
                                            final Set<String> controlSampleSet,
                                            final Set<String> treatSampleSet,
                                            final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                                            final double[] controlMedians,
                                            final double[] treatMedians) {
        final PermutationTTest.Alternative alternative;
        if (svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
            alternative = PermutationTTest.Alternative.GREATER;
        } else if (svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            alternative = PermutationTTest.Alternative.LESS;
        } else {
            throw new IllegalStateException("Unsupported variant type: " + svtype);
        }
        final PermutationTTest.PermTestResult result = PermutationTTest.test(controlMedians, treatMedians, alternative);

        final int numBins = depthMatrix.getNumBins();
        final double[] slicePvals = new double[numBins];
        for (int i = 0; i < numBins; i++) {
            final double[] controlSlice = depthMatrix.slice(i, controlSampleSet);
            final double[] treatSlice = depthMatrix.slice(i, treatSampleSet);
            final PermutationTTest.PermTestResult sliceResult = PermutationTTest.test(controlSlice, treatSlice, alternative);
            slicePvals[i] = sliceResult.pValue();
        }
        Arrays.sort(slicePvals);
        final double secondMaxP;
        if (slicePvals.length > 0) {
            secondMaxP = slicePvals[Math.max(slicePvals.length - 2, 0)];
        } else {
            secondMaxP = Double.NaN;
        }
        return new DepthStats(result.pValue(), secondMaxP);
    }

    private static double combinePValues(final double[] x) {
        if (x.length == 0) {
            return Double.NaN;
        } else if (x.length == 1) {
            return Math.exp(x[0]);
        } else {
            // Fisher's method
            final double chiSq = -2 * MathUtils.sum(x);
            return 1 - new ChiSquaredDistribution(2 * x.length).cumulativeProbability(chiSq);
        }
    }

    protected static double getMedianSeparation(final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                                                final double[] controlMedians, final double[] treatMedians) {
        final int medianSeparationSign;
        if (svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
            medianSeparationSign = 1;
        } else if (svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            medianSeparationSign = -1;
        } else {
            throw new IllegalStateException("Unsupported variant type: " + svtype);
        }
        return medianSeparationSign * (MEDIAN.evaluate(controlMedians) - MEDIAN.evaluate(treatMedians));
    }
}
