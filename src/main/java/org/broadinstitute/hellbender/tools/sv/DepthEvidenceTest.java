package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Genotype;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.jetbrains.annotations.Nullable;
import picard.util.MathUtil;

import java.util.*;
import java.util.stream.Collectors;

public class DepthEvidenceTest {

    protected record SampleSetResult(Set<String> treatSampleSet, Set<String> controlSampleSet) {
    }

    protected record DepthStats(double pValue, double secondMaxP) {
    }

    public record DepthTestResult(double pValue, double secondMaxP, double medianSeparation) {
    }

    protected static final Median MEDIAN = new Median();
    protected final double powerThreshold;

    public DepthEvidenceTest(final double powerThreshold) {
        this.powerThreshold = powerThreshold;
    }

    public DepthTestResult test(final SVCallRecord record, final DepthMatrix depthMatrix) {
        final SampleSetResult sampleSets = getSampleSets(record);
        if (sampleSets == null) return null;
        final double[] controlMedians = getSampleMedians(depthMatrix, sampleSets.controlSampleSet);
        final double[] treatMedians = getSampleMedians(depthMatrix, sampleSets.treatSampleSet);
        final double power = TTestPowerCalculator.powerCalc(controlMedians, treatMedians);
        final GATKSVVCFConstants.StructuralVariantAnnotationType svtype = record.getType();
        final DepthStats result;
        if (!Double.isNaN(power) && power > powerThreshold) {
            result = highPowerTest(depthMatrix, sampleSets.controlSampleSet, sampleSets.treatSampleSet, svtype, controlMedians, treatMedians);
        } else {
            result = lowPowerTest(depthMatrix, sampleSets.controlSampleSet, sampleSets.treatSampleSet, svtype, treatMedians, controlMedians);
        }
        final double medianSeparation = getMedianSeparation(svtype, controlMedians, treatMedians);
        return new DepthTestResult(result.pValue, result.secondMaxP, medianSeparation);
    }

    protected @Nullable SampleSetResult getSampleSets(final SVCallRecord record) {
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
        final Set<String> rawRefSampleSet = Sets.difference(record.getAllSamples(), rawCarrierSampleSet); // TODO: does not account for no-call samples
        final Set<String> eligibleGenotypes = getSamplesWithExpectedCopyNumber(maxExpectedCopyNumber, record.getGenotypes()); // TODO: male-only variants only tested against other males but females could be incorporated here
        final Set<String> treatSampleSet = Sets.intersection(rawCarrierSampleSet, eligibleGenotypes);
        final Set<String> controlSampleSet = Sets.intersection(rawRefSampleSet, eligibleGenotypes);
        return new SampleSetResult(treatSampleSet, controlSampleSet);
    }

    private static int getMaxExpectedCopyNumber(final Collection<Genotype> genotypes) {
        Utils.nonNull(genotypes);
        Utils.nonNull(genotypes);
        final OptionalInt result = genotypes.stream()
                .mapToInt(g -> VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 0))
                .max();
        return result.getAsInt();
    }

    private static Set<String> getSamplesWithExpectedCopyNumber(final int state, final List<Genotype> genotypes) {
        Utils.nonNull(genotypes);
        return genotypes.stream()
                .filter(g -> VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 0) == state)
                .map(Genotype::getSampleName)
                .collect(Collectors.toUnmodifiableSet());
    }

    private double[] getSampleMedians(final DepthMatrix depthMatrix, final Set<String> samples) {
        final double[] medians = new double[samples.size()];
        int i = 0;
        for (final String sample : samples) {
            medians[i++] = MEDIAN.evaluate(depthMatrix.getSample(sample));
        }
        return medians;
    }

    protected DepthStats lowPowerTest(final DepthMatrix depthMatrix,
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
        final NormalDistribution normal = new NormalDistribution();
        final double[] pValues = new double[treatMedians.length];
        final double[] secondMaxPValues = new double[treatMedians.length];
        final double controlMean = MathUtil.mean(controlMedians);
        final double controlStd = MathUtil.stddev(controlMedians, controlMean) * Math.sqrt(controlMedians.length / (double) (controlMedians.length - 1));
        final String[] treatSampleList = treatSampleSet.toArray(new String[0]);
        final int numBins = depthMatrix.getNumBins();
        for (int i = 0; i < treatMedians.length; i++) {
            pValues[i] = Math.log(Math.max(normal.cumulativeProbability(alternativeSign * (treatMedians[i] - controlMean) / controlStd), Double.MIN_VALUE));
            final double[] slicePvals = new double[numBins];
            for (int j = 0; j < numBins; j++) {
                final double[] controlSlice = depthMatrix.slice(j,  controlSampleSet);
                final double[] treatSlice = depthMatrix.slice(j, Collections.singleton(treatSampleList[i]));
                final double controlMeanSlice = MathUtil.mean(controlSlice);
                final double controlStdSlice = MathUtil.stddev(controlSlice, controlMeanSlice) * Math.sqrt(controlSlice.length / (double) (controlSlice.length - 1));
                slicePvals[j] = Math.log(Math.max(normal.cumulativeProbability(alternativeSign * (treatSlice[0] - controlMeanSlice) / controlStdSlice), Double.MIN_VALUE));
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

    protected DepthStats highPowerTest(final DepthMatrix depthMatrix,
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

    private double combinePValues(final double[] x) {
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

    protected double getMedianSeparation(final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                                         final double[] controlMedians, final double[] treatMedians) {
        final int medianSeparationSign;
        if (svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
            medianSeparationSign = 1;
        } else if (svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
            medianSeparationSign = -1;
        } else {
            throw new IllegalStateException("Unsupported variant type: " + svtype);
        }
        final double medianSeparation = medianSeparationSign * (MEDIAN.evaluate(controlMedians) - MEDIAN.evaluate(treatMedians));
        return medianSeparation;
    }
}
