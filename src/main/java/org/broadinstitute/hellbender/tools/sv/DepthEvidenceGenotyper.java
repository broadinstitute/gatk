package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.aggregation.DepthEvidenceTest;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.PloidyTable;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import picard.util.MathUtil;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class DepthEvidenceGenotyper {

    private static final double DEFAULT_COPY_STATE_INCREMENT = 0.5;
    private static final double DEFAULT_COPY_STATE_STD = 0.15; // arbitrary, gives reasonable GQs when training
    private static final int DEFAULT_NUM_STATES = 6; // must be at least 3
    private static final int NO_CARRIER_VARIANT_QUAL = 0;

    private static final NormalDistribution Z_DISTRIBUTION = new NormalDistribution();
    private static final Median MEDIAN = new Median();

    private final SAMSequenceDictionary dictionary;
    private final CopyStateStats[] copyStateStats;
    private final List<String> samples;

    public DepthEvidenceGenotyper(final List<CopyStateStats> cutoffs, final List<String> samples, final SAMSequenceDictionary dictionary) {
        this.dictionary = Utils.nonNull(dictionary);
        this.samples = Utils.nonNull(samples);
        if (cutoffs != null) {
            copyStateStats = cutoffs.toArray(CopyStateStats[]::new);
        } else {
            copyStateStats = new CopyStateStats[DEFAULT_NUM_STATES];
            for (int i = 0; i < copyStateStats.length; i++) {
                copyStateStats[i] = new CopyStateStats(i, i * DEFAULT_COPY_STATE_INCREMENT, DEFAULT_COPY_STATE_STD, (i + 0.5) * DEFAULT_COPY_STATE_INCREMENT);
            }
        }
    }

    private int getCopyState(final double val) {
        for (int i = 0; i < copyStateStats.length; i++) {
            if (val < copyStateStats[i].upperBound) {
                return i;
            }
        }
        return copyStateStats.length - 1;
    }

    public SVCallRecord applyToRecord(final SVCallRecord record,
                                      final DepthGenotypeResult result,
                                      final PloidyTable ploidyTable,
                                      final int maxQuality) {
        Utils.nonNull(record);
        Utils.nonNull(result);
        final Map<String, Object> refinedAttr = new HashMap<>(record.getAttributes());
        refinedAttr.put(GATKSVVCFConstants.DEPTH_VARIANT_QUALITY_ATTRIBUTE, Math.min(result.variantQual, maxQuality));

        final List<Genotype> genotypes = record.getGenotypes();
        final GenotypesContext newGenotypes = GenotypesContext.create(genotypes.size());
        final Map<String, Integer> sampleIndexMap = new HashMap<>();
        for (int i = 0; i < samples.size(); i++) {
            sampleIndexMap.put(samples.get(i), i);
        }
        for (final String sample : samples) {
            final GenotypeBuilder builder;
            final int ploidy;
            if (record.getGenotypes().containsSample(sample)) {
                final Genotype genotype = record.getGenotypes().get(sample);
                builder = new GenotypeBuilder(genotype);
                if (!genotype.hasAnyAttribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT)) {
                    throw new IllegalArgumentException("Sample " + sample + " has no ECN field");
                }
                ploidy = VariantContextGetters.getAttributeAsInt(genotype, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 0);
            } else {
                builder = new GenotypeBuilder(sample);
                if (ploidyTable != null) {
                    ploidy = ploidyTable.get(sample, record.getContigA());
                    builder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, ploidy);
                } else {
                    throw new IllegalArgumentException("Ploidy table required here");
                }
            }
            if (!sampleIndexMap.containsKey(sample)) {
                newGenotypes.add(builder.make());
                continue;
            }
            final DepthGenotype depthGenotype = result.genotypeQuals.get(sampleIndexMap.get(sample));
            if (depthGenotype == null) {
                newGenotypes.add(builder.make());
                continue;
            }
            if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
                builder.alleles(Collections.nCopies(ploidy, Allele.NO_CALL));
            } else if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
                builder.alleles(CanonicalSVCollapser.getDuplicationAllelesFromCopyNumber(record.getRefAllele(), ploidy, depthGenotype.copyState));
            } else if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
                builder.alleles(CanonicalSVCollapser.getDeletionAllelesFromCopyNumber(record.getRefAllele(), ploidy, depthGenotype.copyState));
            } else {
                throw new GATKException("Unsupported variant type here: " + record.getType());
            }
            builder.attribute(GATKSVVCFConstants.DEPTH_COPY_STATE_ATTRIBUTE, depthGenotype.copyState);
            builder.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_QUALITY_ATTRIBUTE, Math.min(Math.round(depthGenotype.quality), maxQuality));
            newGenotypes.add(builder.make());
        }

        // Create new record
        return new SVCallRecord(record.getId(), record.getContigA(), record.getPositionA(),
                record.getStrandA(), record.getContigB(), record.getPositionB(), record.getStrandB(),
                record.getType(), record.getComplexSubtype(), record.getComplexEventIntervals(), record.getLength(),
                record.getEvidence(), record.getAlgorithms(), record.getAlleles(),
                newGenotypes, refinedAttr, record.getFilters(), record.getLog10PError(), dictionary);
    }

    public DepthGenotypeResult genotype(final DepthMatrix depthMatrix) {
        // Get median coverage for each sample
        final double[] medians = DepthEvidenceTest.getSampleMedians(depthMatrix, samples);
        final int numSamples = medians.length;
        if (numSamples == 0) {
            return null;
        }
        // Get genotypes and genotype qualities
        final int[] copyStates = new int[numSamples];
        final List<Double> carrierPValues = new ArrayList<>();
        final List<DepthGenotype> genotypeQuals = new ArrayList<>(numSamples);
        for (int i = 0; i < numSamples; i++) {
            copyStates[i] =  getCopyState(medians[i]);
            final CopyStateStats state = copyStateStats[copyStates[i]];
            final CopyStateStats lowState = copyStateStats[Math.max(copyStates[i] - 1, 0)];
            final CopyStateStats highState = copyStateStats[Math.min(copyStates[i] + 1, copyStateStats.length - 1)];
            final double pValue = getPValue(medians[i], state);
            final double pValueLow = getPValue(medians[i], lowState);
            final double pValueHigh = getPValue(medians[i], highState);
            final double qualLow = getQuality(pValueLow, pValue);
            final double qualHigh = getQuality(pValueHigh, pValue);
            final double qual = copyStates[i] == 0 ? qualHigh : Math.min(qualLow, qualHigh);
            genotypeQuals.add(new DepthGenotype(copyStates[i], qual));
            if (copyStates[i] != 2) {
                // TODO we assume diploid, but this is only used for variant quality
                final double pValueRef = getPValue(medians[i], copyStateStats[2]);
                carrierPValues.add(pValueRef);
            }
        }
        // Get variant quality
        final double variantQual;
        if (carrierPValues.isEmpty()) {
            variantQual = NO_CARRIER_VARIANT_QUAL;
        } else {
            final double median = MEDIAN.evaluate(carrierPValues.stream().mapToDouble(Double::doubleValue).toArray());
            variantQual = -10 * Math.log10(median);
        }
        return new DepthGenotypeResult(medians, copyStates, genotypeQuals, variantQual);
    }

    public List<CopyStateStats> train(final Collection<DepthGenotypeResult> genotypes, final int numTrainingStates) {
        final int[] copyStates = genotypes.stream().map(DepthGenotypeResult::copyStates).flatMapToInt(Arrays::stream).toArray();
        final double[] sampleDepths = genotypes.stream().map(DepthGenotypeResult::sampleDepths).flatMapToDouble(Arrays::stream).toArray();
        final List<CopyStateStats> trained = new ArrayList<>(numTrainingStates);
        final double[] means = new double[numTrainingStates];
        final double[] stdDevs = new double[numTrainingStates];
        for (int i = 0; i < numTrainingStates; i++) {
            final List<Integer> samples = new ArrayList<>();
            for (int j = 0; j < sampleDepths.length; j++) {
                if (copyStates[j] == i) {
                    samples.add(j);
                }
            }
            if (samples.isEmpty()) {
                means[i] = DEFAULT_COPY_STATE_INCREMENT * i;
                stdDevs[i] = i == 0 ? DEFAULT_COPY_STATE_INCREMENT * 0.5 : stdDevs[i - 1];
            } else {
                final double[] depths = new double[samples.size()];
                for (int j = 0; j < samples.size(); j++) {
                    depths[j] = sampleDepths[samples.get(j)];
                }
                means[i] = MathUtil.mean(depths);
                stdDevs[i] = MathUtil.stddev(depths, means[i]);
            }
        }
        for (int i = 0; i < numTrainingStates; i++) {
            final double upperBound;
            if (i == numTrainingStates - 1) {
                upperBound = (i + 0.5) * DEFAULT_COPY_STATE_INCREMENT;
            } else {
                upperBound = getCutoff(means[i], stdDevs[i], means[i + 1], stdDevs[i + 1]);
            }
            trained.add(new CopyStateStats(i, means[i], stdDevs[i], upperBound));
        }
        return trained;
    }

    private static double getCutoff(final double mean1, final double std1, final double mean2, final double std2) {
        final double w1 = 1. / std1;
        final double w2 = 1. / std2;
        return (w1 * mean1 + w2 * mean2) / (w1 + w2);
    }

    private static double getQuality(final double p1, final double p2) {
        if (Double.isNaN(p1) || Double.isNaN(p2)) {
            return Double.NaN;
        }
        Utils.validateArg(0 <= p1 && 1 >= p1, "p1 not on [0, 1]");
        Utils.validateArg(0 <= p2 && 1 >= p2, "p2 not on [0, 1]");
        Utils.validateArg(p1 <= p2, "p1 cannot be greater than p2");
        if (p2 == 0) {
            // p1 also 0
            return 0;
        } else if (p1 == 0) {
            return Double.POSITIVE_INFINITY;
        }
        return -10. * Math.log10(p1 / p2);
    }

    private static double getPValue(final double median, final CopyStateStats state) {
        return 1. - Z_DISTRIBUTION.cumulativeProbability(Math.abs(state.mean - median) / state.stdDev);
    }

    public record CopyStateStats(int copyState, double mean, double stdDev, double upperBound) {}
    public record DepthGenotypeResult(double[] sampleDepths, int[] copyStates, List<DepthGenotype> genotypeQuals, double variantQual) {}
    public record DepthGenotype(int copyState, double quality) {}
}
