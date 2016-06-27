package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * Genotype predetermined segments passed in the inputs together with the targets and their coverage per sample.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        programGroup = CopyNumberProgramGroup.class,
        summary = "Genotype locations for copy number variation in germline samples using a HMM",
        oneLineSummary = "Genotype location for copy number variation"
)
public final class GenotypeCopyNumberTriStateSegments extends CopyNumberTriStateSegmentCaller {

    public static final String DISCOVERY_KEY = "DSCVR";
    public static final String NUMBER_OF_TARGETS_KEY = "NTARGETS";
    public static final String SOME_QUALITY_KEY = "SQ";
    public static final String START_QUALITY_KEY = "LQ";
    public static final String END_QUALITY_KEY = "RQ";

    public static final String DISCOVERY_FILE_SHORT_NAME = "segments";
    public static final String DISCOVERY_FILE_FULL_NAME = "segmentsFile";

    public static final String DISCOVERY_TRUE = "Y";
    public static final String DISCOVERY_FALSE = "N";

    public static final int MAX_GQ = 99;

    @Argument(
            doc = "Discovered segments file",
            fullName = DISCOVERY_FILE_FULL_NAME,
            shortName = DISCOVERY_FILE_SHORT_NAME
    )
    protected File segmentsFile;

    private VariantContextWriter outputWriter;

    @Override
    protected void openOutput(final File outputFile, final CopyNumberTriStateHiddenMarkovModel model,
                              final TargetCollection<Target> targets, ReadCountCollection inputCounts) {
        outputWriter = GATKVariantContextUtils.createVCFWriter(outputFile, null, false);
        outputWriter.writeHeader(composeHeader(inputCounts.columnNames()));
    }

    private VCFHeader composeHeader(final List<String> sampleNames) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), sampleNames);

        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(), VCFHeaderVersion.VCF4_2.getVersionString()));
        CopyNumberTriStateAllele.addHeaderLinesTo(result);
        result.addMetaDataLine(new VCFHeaderLine("command",this.getCommandLine()));

        // FORMAT.
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.Integer, "Genotype"));
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_PL_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"));
        result.addMetaDataLine(new VCFFormatHeaderLine(DISCOVERY_KEY, 1, VCFHeaderLineType.Character, "Y if this segment was discovered in that sample, N otherwise"));
        result.addMetaDataLine(new VCFFormatHeaderLine(SOME_QUALITY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Quality of the region to contain at least one target with that call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(START_QUALITY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Quality of the segment to start exactly at that location for alternative allele calls only"));
        result.addMetaDataLine(new VCFFormatHeaderLine(END_QUALITY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Quality of the segment to end exactly at that location for alternative allele calls only"));
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1, VCFHeaderLineType.Integer, "Genotype call quality as the difference between the best and second best PL"));
        // INFO
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.ALLELE_COUNT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Allele count in genotypes, for each ALT allele, in the same order as listed"));
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele Frequency, for each ALT allele, in the same order as listed"));
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.ALLELE_NUMBER_KEY, 1, VCFHeaderLineType.Integer, "Total number of alleles in called genotypes"));
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1, VCFHeaderLineType.Integer, "End coordinate of this variant"));
        result.addMetaDataLine(new VCFInfoHeaderLine(NUMBER_OF_TARGETS_KEY, 1, VCFHeaderLineType.Integer, "Number of targets enclosed in this variant coordinates"));
        return result;
    }

    @Override
    protected void closeOutput(final File outputFile) {
        outputWriter.close();
    }

    @Override
    protected void makeCalls(final CopyNumberTriStateHiddenMarkovModel model, final TargetCollection<Target> targets, final ReadCountCollection inputCounts) {
        logger.info("Composing list of segment intervals to genotype ...");
        final List<GenotypingSegment> segments = composeGenotypingSegments(segmentsFile, targets);
        logger.info(String.format("A total of %d segments to genotype found", segments.size()));
        final List<ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState>> fbResults = runFWBWAlgorithm(model, inputCounts);
        for (final GenotypingSegment segment : segments) {
            final VariantContext variant = composeVariantContext(segment, inputCounts.columnNames(), fbResults);
            outputWriter.add(variant);
        }
    }


    private VariantContext composeVariantContext(final GenotypingSegment segment,
                                                 final List<String> sampleNames,
                                                 final List<ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState>> fbResults) {
        final VariantContextBuilder builder = new VariantContextBuilder();
        builder.alleles(CopyNumberTriStateAllele.PLAIN_ALL_ALLELES);
        builder.chr(segment.getContig());
        builder.start(segment.getStart());
        builder.stop(segment.getEnd());
        builder.id(String.format("CNV_%s_%d_%d", segment.getContig(), segment.getStart(), segment.getEnd()));
        final IndexRange targetIndexes = segment.getTargetIndexes();
        final List<Genotype> genotypes = IntStream.range(0, sampleNames.size())
                .mapToObj(i -> {
                    final String sample = sampleNames.get(i);
                    final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult = fbResults.get(i);
                    final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
                    final double[] log10GP = calculateLog10GP(targetIndexes, fbResult);
                    final GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(log10GP);
                    final double[] SQ = calculateSQ(targetIndexes, fbResult);
                    final double[] LQ = calculateLQ(targetIndexes, fbResult);
                    final double[] RQ = calculateRQ(targetIndexes, fbResult);
                    final int[] PL = likelihoods.getAsPLs();
                    final int GQ = calculateGQ(PL);
                    final int genotypeCall = MathUtils.maxElementIndex(log10GP);
                    genotypeBuilder.PL(PL);
                    genotypeBuilder.alleles(Collections.singletonList(CopyNumberTriStateAllele.ALL_ALLELES.get(genotypeCall)));
                    genotypeBuilder.attribute(DISCOVERY_KEY, segment.containingSamples.contains(sample) ? DISCOVERY_TRUE : DISCOVERY_FALSE);
                    genotypeBuilder.attribute(SOME_QUALITY_KEY, SQ);
                    genotypeBuilder.attribute(START_QUALITY_KEY, LQ);
                    genotypeBuilder.attribute(END_QUALITY_KEY, RQ);
                    genotypeBuilder.GQ(GQ);
                    return genotypeBuilder.make();
                }).collect(Collectors.toList());
        final int alleleNumber = genotypes.size();
        final int deletionCount = (int) genotypes.stream()
                .filter(g -> g.getAllele(0) == CopyNumberTriStateAllele.DEL).count();
        final int duplicationCount = (int) genotypes.stream()
                .filter(g -> g.getAllele(0) == CopyNumberTriStateAllele.DUP).count();
        final double deletionFrequency = deletionCount / (double) alleleNumber;
        final double duplicationFrequency = duplicationCount / (double) alleleNumber;

        builder.attribute(VCFConstants.END_KEY, segment.getEnd());
        builder.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] { deletionFrequency, duplicationFrequency});
        builder.attribute(VCFConstants.ALLELE_COUNT_KEY, new int[] { deletionCount, duplicationCount});
        builder.attribute(VCFConstants.ALLELE_NUMBER_KEY, alleleNumber);
        builder.attribute(NUMBER_OF_TARGETS_KEY, targetIndexes.size());
        builder.genotypes(genotypes);
        return builder.make();
    }

    private int calculateGQ(final int[] PL) {
        int best = PL[0];
        int secondBest = Integer.MAX_VALUE;
        for (int i = 1; i < PL.length; i++) {
            final int value = PL[i];
            if (value <= best) {
                secondBest = best;
                best = value;
            } else if (value < secondBest) {
                secondBest = value;
            }
        }
        return Math.min(secondBest - best, MAX_GQ);
    }

    private double[] calculateRQ(final IndexRange targetIndexes, final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        return CopyNumberTriStateAllele.ALTERNATIVE_ALLELES.stream()
                .mapToDouble(allele ->
                        logProbToPhredScore(logEndProbability(targetIndexes.to - 1,
                                allele.state, fbResult), true)).toArray();
    }

    private double[] calculateLQ(final IndexRange targetIndexes, final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        return CopyNumberTriStateAllele.ALTERNATIVE_ALLELES.stream()
                .mapToDouble(allele ->
                    logProbToPhredScore(logStartProbability(targetIndexes.from,
                            allele.state, fbResult), true)).toArray();
    }

    private double[] calculateSQ(final IndexRange targetIndexes, final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        return CopyNumberTriStateAllele.ALTERNATIVE_ALLELES.stream()
            .mapToDouble(allele ->
                    logProbToPhredScore(logSomeProbability(targetIndexes.from, targetIndexes.size(),
                            allele.state, fbResult), true)).toArray();
    }

    private double[] calculateLog10GP(final IndexRange targetIndexes, final ForwardBackwardAlgorithm.Result<Double, Target, CopyNumberTriState> fbResult) {
        return CopyNumberTriStateAllele.ALL_ALLELES.stream()
            .mapToDouble(allele -> fbResult.logProbability(targetIndexes.from, targetIndexes.to, allele.state) * INV_LN_10)
            .map(CopyNumberTriStateSegmentCaller::roundPhred).toArray();
    }

    private List<ForwardBackwardAlgorithm.Result<Double,Target,CopyNumberTriState>> runFWBWAlgorithm(final CopyNumberTriStateHiddenMarkovModel model, final ReadCountCollection inputCounts) {
        logger.info(String.format("Running the forward-backward algorithm in a total of %d samples over %d targets", inputCounts.columnNames().size(), inputCounts.targets().size()));
        return IntStream.range(0, inputCounts.columnNames().size())
                .mapToObj(i -> {
                    final List<Double> data = DoubleStream.of(inputCounts.counts().getColumn(i)).boxed().collect(Collectors.toList());
                    return ForwardBackwardAlgorithm.apply(data, inputCounts.targets(), model);
                })
                .collect(Collectors.toList());
    }


    private List<GenotypingSegment> composeGenotypingSegments(final File segmentsFile, final TargetCollection<Target> targets) {

        try (final CopyNumberTriStateSegmentRecordReader reader = new CopyNumberTriStateSegmentRecordReader(segmentsFile)) {
            final List<CopyNumberTriStateSegmentRecord> records = reader.stream()
                    .filter(s -> s.getSegment().getCall() != CopyNumberTriState.NEUTRAL)
                    .sorted(Comparator.comparing(
                            CopyNumberTriStateSegmentRecord::getSegment,
                            IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR))
                    .collect(Collectors.toList());
            final List<GenotypingSegment> result = new ArrayList<>();
            for (final CopyNumberTriStateSegmentRecord record : records) {
                if (result.isEmpty()) {
                    result.add(new GenotypingSegment(record, targets));
                } else {
                    final GenotypingSegment lastSegment = result.get(result.size() - 1);
                    if (lastSegment.getInterval().equals(record.getSegment().getInterval())) {
                        lastSegment.addSample(record.getSampleName());
                    } else {
                        result.add(new GenotypingSegment(record, targets));
                    }
                }
            }
            return result;
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(segmentsFile, ex);
        }
    }

    /**
     * Use to collect information about segment to be genotyped.
     */
    final class GenotypingSegment implements Locatable {

        private final SimpleInterval interval;
        private final IndexRange targetIndexes;
        private final Set<String> containingSamples;

        public GenotypingSegment(final CopyNumberTriStateSegmentRecord record, final TargetCollection<Target> targets) {
            this(record.getSegment().getInterval(), targets.indexRange(record.getSegment()));
            addSample(record.getSampleName());
        }

        public GenotypingSegment(final SimpleInterval interval, final IndexRange targetIndexes) {
            this.interval = Utils.nonNull(interval);
            this.targetIndexes = Utils.nonNull(targetIndexes);
            this.containingSamples = new HashSet<>();
        }

        public void addSample(final String name) {
            this.containingSamples.add(name);
        }

        public IndexRange getTargetIndexes() {
            return targetIndexes;
        }

        public SimpleInterval getInterval() {
            return interval;
        }

        public int getTargetCount() {
            return targetIndexes.size();
        }

        @Override
        public String getContig() {
            return interval.getContig();
        }

        @Override
        public int getStart() {
            return interval.getStart();
        }

        @Override
        public int getEnd() {
            return interval.getEnd();
        }
    }
}
