package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledContigPloidyCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntervalCopyNumberGenotypingData;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVIntervalVariantDecoder;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

public class DepthEvidenceAggregator {

    private final List<VCFFileReader> posteriorsReaders;
    private final SAMSequenceDictionary dictionary;
    private final List<String> samples;
    private final GermlineCNVIntervalVariantDecoder cnvDecoder;
    private final List<IntegerCopyNumberState> copyStates;
    private final int numCopyStates;
    private final double copyNeutralPrior;
    private Map<String,Map<String,Integer>> sampleContigPloidyMap;

    private String currentContig;
    private List<IntervalTree<Map<String,double[]>>> currentPosteriorsTreeList;

    public static final int COPY_NEUTRAL_PRIOR_BASIS_LENGTH = 1000;

    public DepthEvidenceAggregator(final List<VCFFileReader> posteriorsReaders,
                                   final Collection<CalledContigPloidyCollection> contigPloidyCollections,
                                   final double copyNeutralPrior,
                                   final List<String> samples,
                                   final SAMSequenceDictionary dictionary) {
        Utils.nonNull(posteriorsReaders);
        Utils.nonNull(contigPloidyCollections);
        Utils.nonNull(samples);
        Utils.nonNull(dictionary);
        this.posteriorsReaders = posteriorsReaders;
        this.samples = samples;
        this.dictionary = dictionary;
        this.copyNeutralPrior = copyNeutralPrior;

        cnvDecoder = new GermlineCNVIntervalVariantDecoder(contigPloidyCollections);
        sampleContigPloidyMap = cnvDecoder.getSampleContigPloidyMap();

        final VariantContext exampleVariant = posteriorsReaders.get(0).iterator().next();
        copyStates = getCopyNumberStates(exampleVariant);
        numCopyStates = copyStates.size();
        validateCopyStates();
    }

    public VariantContext apply(final SVCallRecord call, final VariantContext baseVariant) {
        Utils.nonNull(call);
        Utils.nonNull(baseVariant);
        final Map<String, CopyNumberPosteriorDistribution> variantPosteriors = samples.stream()
                .map(s -> new Tuple2<>(s, getCallPosterior(call, s)))
                .filter(t -> t._2 != null)
                .collect(Collectors.toMap(t -> t._1, t -> t._2));
        return createBaseVariant(baseVariant, variantPosteriors);
    }

    public CopyNumberPosteriorDistribution getCallPosterior(final SVCallRecord call, final String sample) {
        Utils.nonNull(call);
        Utils.nonNull(sample);
        if (!call.getContigA().equals(currentContig)) {
            currentContig = call.getContigA();
            currentPosteriorsTreeList = posteriorsReaders.stream().map(this::getCurrentPosteriorsTree).collect(Collectors.toList());
        }
        if (!(call.getType().equals(StructuralVariantType.DEL)  || call.getType().equals(StructuralVariantType.DUP)
                || call.getType().equals(StructuralVariantType.BND))) {
            return null;
        }
        if (!call.getContigA().equals(call.getContigB())) {
            return null;
        }
        final SimpleInterval interval = new SimpleInterval(call.getContigA(), call.getPositionA(), call.getPositionB());
        final List<IntervalTree.Node<Map<String,double[]>>> posteriorsList = getBestPosteriors(interval);
        return getIntervalPosterior(interval, posteriorsList, sample);
    }

    public VariantContext createBaseVariant(final VariantContext variant,
                                            final Map<String, CopyNumberPosteriorDistribution> variantPosteriors) {
        Utils.nonNull(variant);
        Utils.nonNull(variantPosteriors);
        final VariantContextBuilder variantBuilder = new VariantContextBuilder(variant);
        final Map<String, List<Integer>> copyStateQuals = getCopyStateQuals(variantPosteriors);
        final List<Genotype> newGenotypes = new ArrayList<>(variant.getGenotypes().size());
        for (final Genotype genotype : variant.getGenotypes()) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);
            final String sample = genotype.getSampleName();
            genotypeBuilder.attribute(GATKSVVCFConstants.COPY_NUMBER_LOG_POSTERIORS_KEY, copyStateQuals.get(sample));
            genotypeBuilder.attribute(GATKSVVCFConstants.NEUTRAL_COPY_NUMBER_KEY, getSamplePloidy(sample, variant.getContig()));
            newGenotypes.add(genotypeBuilder.make());
        }
        variantBuilder.genotypes(newGenotypes);
        return variantBuilder.make();
    }

    private Map<String, List<Integer>> getCopyStateQuals(final Map<String, CopyNumberPosteriorDistribution> variantPosteriors) {
        final Map<String, List<Integer>> copyStateQuals = new HashMap<>(SVUtils.hashMapCapacity(variantPosteriors.size()));
        for (final String sample : variantPosteriors.keySet()) {
            final CopyNumberPosteriorDistribution dist = variantPosteriors.get(sample);
            final List<Integer> copyStatePhred = dist.getIntegerCopyNumberStateList().stream()
                    .map(dist::getCopyNumberPosterior)
                    .map(QualityUtils::logProbToPhred)
                    .map(Math::round)
                    .map(Long::intValue)
                    .collect(Collectors.toList());
            copyStateQuals.put(sample, copyStatePhred);
        }
        return copyStateQuals;
    }

    private void validateCopyStates() {
        for (int i = 1; i < posteriorsReaders.size(); i++) {
            final List<IntegerCopyNumberState> otherCopyStates = getCopyNumberStates(posteriorsReaders.get(i).iterator().next());
            if (!copyStates.equals(otherCopyStates)) {
                throw new UserException.BadInput("CNV VCFs do not contain identical copy number states.");
            }
        }
    }

    private List<IntegerCopyNumberState> getCopyNumberStates(final VariantContext variant) {
        final SimpleInterval interval = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getEnd());
        return cnvDecoder.parseGenotype(variant.getGenotype(0), interval)
                .getCopyNumberPosteriorDistribution().getIntegerCopyNumberStateList();
    }

    private IntervalTree<Map<String,double[]>> getCurrentPosteriorsTree(final VCFFileReader reader) {
        final SAMSequenceRecord contigRecord = dictionary.getSequence(currentContig);
        if (contigRecord == null) {
            throw new UserException.MissingContigInSequenceDictionary(currentContig, dictionary);
        }
        final Iterator<VariantContext> posteriorsIter = reader.query(currentContig, 1, contigRecord.getSequenceLength());
        final IntervalTree<Map<String,double[]>> tree = new IntervalTree<>();
        while (posteriorsIter.hasNext()) {
            final VariantContext variant = posteriorsIter.next();
            final SimpleInterval interval = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getEnd());
            final Map<String,double[]> samplePosteriorMap = new HashMap<>(SVUtils.hashMapCapacity(samples.size()));
            for (final Genotype genotype : variant.getGenotypes()) {
                final IntervalCopyNumberGenotypingData data = cnvDecoder.parseGenotype(genotype, interval);
                final double[] posteriors = new double[numCopyStates];
                int i = 0;
                for (final IntegerCopyNumberState state : copyStates) {
                    posteriors[i] = data.getCopyNumberPosteriorDistribution().getCopyNumberPosterior(state);
                    i++;
                }
                samplePosteriorMap.put(genotype.getSampleName(), posteriors);
            }
            tree.put(interval.getStart(), interval.getEnd(), samplePosteriorMap);
        }
        return tree;
    }

    private List<IntervalTree.Node<Map<String,double[]>>> getBestPosteriors(final SimpleInterval interval) {
        int maxOverlap = -1;
        int maxIntervalCount = 0;
        List<IntervalTree.Node<Map<String,double[]>>> maxNodeList = null;
        for (int i = 0; i < currentPosteriorsTreeList.size(); i++) {
            final List<IntervalTree.Node<Map<String,double[]>>> nodeList = Lists.newArrayList(currentPosteriorsTreeList.get(i)
                    .overlappers(interval.getStart(), interval.getEnd()));
            final int overlap = nodeList.stream()
                    .map(node -> new SimpleInterval(interval.getContig(), node.getStart(), node.getEnd()))
                    .mapToInt(nodeInterval -> nodeInterval.intersect(interval).getLengthOnReference())
                    .sum();
            if (overlap > maxOverlap || (overlap == maxOverlap && nodeList.size() > maxIntervalCount)) {
                maxOverlap = overlap;
                maxIntervalCount = nodeList.size();
                maxNodeList = nodeList;
            }
        }
        return maxNodeList;
    }

    private CopyNumberPosteriorDistribution getIntervalPosterior(final SimpleInterval variantInterval,
                                                                 final List<IntervalTree.Node<Map<String,double[]>>> posteriorsList,
                                                                 final String sample) {
        final double[] copyStateSums = new double[numCopyStates];
        Arrays.fill(copyStateSums, Double.MIN_VALUE);
        int overlapSize = 0;
        for (final IntervalTree.Node<Map<String,double[]>> node : posteriorsList) {
            final SimpleInterval posteriorInterval = new SimpleInterval(currentContig, node.getStart(), node.getEnd());
            final int overlap = posteriorInterval.intersect(variantInterval).size();
            final double overlapFraction = overlap / (double) posteriorInterval.getLengthOnReference();
            overlapSize += overlap;
            final double[] dist = node.getValue().get(sample);
            for (int j = 0; j < numCopyStates; j++) {
                copyStateSums[j] += dist[j] * overlapFraction;
            }
        }

        // Fill in missing copy number posterior intervals with a prior
        final double unsupportedIntervals = (variantInterval.size() - overlapSize) / (double) COPY_NEUTRAL_PRIOR_BASIS_LENGTH;
        if (unsupportedIntervals > 0) {
            final int ploidy = getSamplePloidy(sample, variantInterval.getContig());
            final double logNeutralProb = FastMath.log(copyNeutralPrior);
            final double logNonNeutralProb = FastMath.log((1.0 - copyNeutralPrior) / (copyStateSums.length - 1));
            for (int i = 0; i < copyStateSums.length; i++) {
                if (i != ploidy) {
                    copyStateSums[i] += logNonNeutralProb * unsupportedIntervals;
                } else {
                    copyStateSums[i] += logNeutralProb * unsupportedIntervals;
                }
            }
        }

        double denom = 0;
        final double maxStateSum = DoubleStream.of(copyStateSums).max().getAsDouble();
        for (int i = 0; i < copyStateSums.length; i++) {
            // Normalize to avoid underflow error
            copyStateSums[i] -= maxStateSum;
            denom += FastMath.exp(copyStateSums[i]);
        }
        final double logDenom = Math.log(denom);
        final Map<IntegerCopyNumberState,Double> eventPosterior = new HashMap<>(SVUtils.hashMapCapacity(numCopyStates));
        for (int i = 0; i < copyStateSums.length; i++) {
            final Double p = copyStateSums[i] - logDenom;
            eventPosterior.put(new IntegerCopyNumberState(i), p);
        }
        return new CopyNumberPosteriorDistribution(eventPosterior);
    }

    private int getSamplePloidy(final String sample, final String contig) {
        if (!sampleContigPloidyMap.containsKey(sample)) {
            throw new UserException("Could not find ploidy calls for sample: " + sample);
        }
        final Map<String,Integer> contigPloidyMap = sampleContigPloidyMap.get(sample);
        if (!contigPloidyMap.containsKey(contig)) {
            throw new UserException("Could not find ploidy calls for contig: " + contig);
        }
        return contigPloidyMap.get(contig);
    }
}
