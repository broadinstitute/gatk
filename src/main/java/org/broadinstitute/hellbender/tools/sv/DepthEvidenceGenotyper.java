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
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import picard.util.MathUtil;

import java.util.*;
import java.util.function.Function;

public class DepthEvidenceGenotyper {

    private static final double DEFAULT_COPY_STATE_INCREMENT = 0.5;
    private static final double DEFAULT_CUTOFF_OFFSET = DEFAULT_COPY_STATE_INCREMENT / 2;
    private static final double DEFAULT_COPY_STATE_STD = 0.15; // arbitrary, gives reasonable GQs when training
    private static final int DEFAULT_NUM_STATES = 6; // must be at least 3
    private static final int NO_CARRIER_VARIANT_QUAL = 0;
    private static final int NO_DATA_VARIANT_QUAL = 0;
    private static final int NO_DATA_COPY_STATE = 2;
    private static final int NO_DATA_GENOTYPE_QUAL = 2;

    private static final NormalDistribution Z_DISTRIBUTION = new NormalDistribution();
    private static final Median MEDIAN = new Median();

    private final SAMSequenceDictionary dictionary;
    private final CopyStateStats[] copyStateStats;
    private final List<String> samples;
    private final int maxQual;

    public DepthEvidenceGenotyper(final List<CopyStateStats> cutoffs, final List<String> samples, final int maxQual, final SAMSequenceDictionary dictionary) {
        this.dictionary = Utils.nonNull(dictionary);
        this.samples = Utils.nonNull(samples);
        this.maxQual = maxQual;
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
        return (int) ((val + DEFAULT_CUTOFF_OFFSET) / DEFAULT_COPY_STATE_INCREMENT);
    }

    public List<String> getSamplesInOrder() {
        return samples;
    }

    /*
    public SVCallRecord applyToRecord(final SVCallRecord record,
                                      final DepthGenotypeResult result,
                                      final int maxQuality) {
        Utils.nonNull(record);
        Utils.nonNull(result);
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
                throw new IllegalArgumentException("Sample " + sample + " is in the samples list but does not have a genotype field");
            }
            if (!sampleIndexMap.containsKey(sample)) {
                newGenotypes.add(builder.make());
                continue;
            }
            final int copyState = result.copyStates[sampleIndexMap.get(sample)];
            final double quality = result.genotypeQuals[sampleIndexMap.get(sample)];
            if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
                builder.alleles(Collections.nCopies(ploidy, Allele.NO_CALL));
            } else if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DUP) {
                builder.alleles(CanonicalSVCollapser.getDuplicationAllelesFromCopyNumber(record.getRefAllele(), ploidy, copyState));
            } else if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.DEL) {
                builder.alleles(CanonicalSVCollapser.getDeletionAllelesFromCopyNumber(record.getRefAllele(), ploidy, copyState));
            } else {
                throw new GATKException("Unsupported variant type here: " + record.getType());
            }
            builder.attribute(GATKSVVCFConstants.DEPTH_COPY_STATE_ATTRIBUTE, copyState);
            builder.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_QUALITY_ATTRIBUTE, Math.min(Math.round(quality), maxQuality));
            newGenotypes.add(builder.make());
        }

        // Create new record
        return new SVCallRecord(record.getId(), record.getContigA(), record.getPositionA(),
                record.getStrandA(), record.getContigB(), record.getPositionB(), record.getStrandB(),
                record.getType(), record.getComplexSubtype(), record.getComplexEventIntervals(), record.getLength(),
                record.getEvidence(), record.getAlgorithms(), record.getAlleles(),
                newGenotypes, record.getAttributes(), record.getFilters(), record.getLog10PError(), dictionary);
    }
     */

    private CopyStateStats getCopyStateStats(final int copyState) {
        if (copyState < copyStateStats.length) {
            return copyStateStats[copyState];
        }
        final double mean = DEFAULT_COPY_STATE_INCREMENT * copyState;
        return new CopyStateStats(copyState, mean, DEFAULT_COPY_STATE_STD, mean + DEFAULT_CUTOFF_OFFSET);
    }

    public DepthGenotypeResult genotype(final DepthMatrix depthMatrix) {
        if (depthMatrix.getNumBins() == 0) {
            // Edge case where there are no overlapping read depth bins
            final double[] medians = new double[samples.size()];
            Arrays.fill(medians, Double.NaN);
            final int[] copyStates = new int[samples.size()];
            Arrays.fill(copyStates, NO_DATA_COPY_STATE);
            final double[] qualities = new double[samples.size()];
            Arrays.fill(qualities, NO_DATA_GENOTYPE_QUAL);
            return new DepthGenotypeResult(medians, copyStates, qualities, NO_DATA_VARIANT_QUAL);
        }
        // Get median coverage for each sample
        final double[] medians = DepthEvidenceTest.getSampleMedians(depthMatrix, samples);
        final int numSamples = medians.length;
        if (numSamples == 0) {
            return null;
        }
        // Get genotypes and genotype qualities
        final int[] copyStates = new int[numSamples];
        final List<Double> carrierPValues = new ArrayList<>();
        final double[] genotypeQuals = new double[numSamples];
        // Reference copy state, whose std dev is used for GQ calcs
        final CopyStateStats refCopyState = getCopyStateStats(2);
        for (int i = 0; i < numSamples; i++) {
            copyStates[i] =  getCopyState(medians[i]);
            final CopyStateStats state = getCopyStateStats(copyStates[i]);
            final CopyStateStats lowState = getCopyStateStats(Math.max(copyStates[i] - 1, 0));
            final CopyStateStats highState = getCopyStateStats(copyStates[i] + 1);
            final double pValue = getPValue(medians[i], state, refCopyState);
            final double pValueLow = getPValue(medians[i], lowState, refCopyState);
            final double pValueHigh = getPValue(medians[i], highState, refCopyState);
            final double qualLow = getQuality(pValueLow, pValue, maxQual);
            final double qualHigh = getQuality(pValueHigh, pValue, maxQual);
            genotypeQuals[i] = copyStates[i] == 0 ? qualHigh : Math.min(qualLow, qualHigh);
            if (copyStates[i] != 2) {
                // TODO we assume diploid, but this is only used for variant quality
                final double pValueRef = getPValue(medians[i], refCopyState, refCopyState);
                carrierPValues.add(pValueRef);
            }
        }
        // Get variant quality
        final double variantQual;
        if (carrierPValues.isEmpty()) {
            variantQual = NO_CARRIER_VARIANT_QUAL;
        } else {
            final double median = MEDIAN.evaluate(carrierPValues.stream().mapToDouble(Double::doubleValue).toArray());
            variantQual = Math.min(-10 * Math.log10(median), maxQual);
        }
        // TODO : implement --geno_adjust option for shifted variants
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

    private static double getQuality(final double p1, final double p2, final double max) {
        if (Double.isNaN(p1) || Double.isNaN(p2)) {
            return Double.NaN;
        }
        Utils.validateArg(0 <= p1 && 1 >= p1, "p1 not on [0, 1]");
        Utils.validateArg(0 <= p2 && 1 >= p2, "p2 not on [0, 1]");
        if (p1 > p2) {
            // case where the assigned copy state contradicts likelihoods
            return 0;
        }
        if (p2 == 0) {
            // p1 also 0
            return 0;
        } else if (p1 == 0) {
            return max;
        }
        return Math.min(-10. * Math.log10(p1 / p2), max);
    }

    private static double getPValue(final double median, final CopyStateStats state, final CopyStateStats refState) {
        return 1. - Z_DISTRIBUTION.cumulativeProbability(Math.abs(state.mean - median) / refState.stdDev);
    }

    public record CopyStateStats(int copyState, double mean, double stdDev, double upperBound) {}
    public record DepthGenotypeResult(double[] sampleDepths, int[] copyStates, double[] genotypeQuals, double variantQual) {}

    public static final class DepthTableParser {

        private static final String COPY_STATE_COLUMN = "copy_state";
        private static final String MEAN_COLUMN = "mean";
        private static final String STD_DEV_COLUMN = "sd";
        private static final String CUTOFFS_COLUMN = "cutoffs";
        public static final TableColumnCollection CUTOFFS_COLUMNS = new TableColumnCollection(Arrays.asList(COPY_STATE_COLUMN, MEAN_COLUMN, STD_DEV_COLUMN, CUTOFFS_COLUMN));

        public void composeCutoffsLine(DepthEvidenceGenotyper.CopyStateStats stats, DataLine dataLine) {
            dataLine.append(stats.copyState());
            dataLine.append(stats.mean());
            dataLine.append(stats.stdDev());
            dataLine.append(stats.upperBound());
        }

        public Function<DataLine, CopyStateStats> tableParser(TableColumnCollection columns, Function<String, RuntimeException> exceptionFactory) {
            // Check for expected columns
            for (final String column : CUTOFFS_COLUMNS.names()) {
                if (!columns.contains(column)) {
                    throw exceptionFactory.apply("Missing column " + column);
                }
            }
            // Check there are no extra columns
            if (columns.columnCount() != CUTOFFS_COLUMNS.columnCount()) {
                throw exceptionFactory.apply("Expected " + columns.columnCount() + " columns but found " + columns.columnCount());
            }
            return this::parseTableLine;
        }

        public DepthEvidenceGenotyper.CopyStateStats parseTableLine(final DataLine dataLine) {
            final int copyState = Integer.parseInt(dataLine.get(COPY_STATE_COLUMN));
            final double mean = Double.parseDouble(dataLine.get(MEAN_COLUMN));
            final double stdDev = Double.parseDouble(dataLine.get(STD_DEV_COLUMN));
            final double cutoff = Double.parseDouble(dataLine.get(CUTOFFS_COLUMN));
            return new DepthEvidenceGenotyper.CopyStateStats(copyState, mean, stdDev, cutoff);
        }

    }
}
