package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledContigPloidyCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.LocatableCopyNumberPosteriorDistributionCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntervalCopyNumberGenotypingData;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.LocatableCopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVIntervalVariantDecoder;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

/**
 * Calculates copy number posteriors for a given set of structural variants. Supports multiple samples.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Structural variant VCF
 *     </li>
 *     <li>
 *         Germline copy number intervals VCF
 *     </li>
 *     <li>
 *         Germline contig ploidy calls
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Structural variant VCF with copy number posteriors
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVCopyNumberPosteriors
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Collects read counts at specified intervals",
        oneLineSummary = "Collects read counts at specified intervals",
        programGroup = CoverageAnalysisProgramGroup.class
)
@DocumentedFeature
public final class SVCopyNumberPosteriors extends VariantWalker {
    public static final String COPY_NUMBER_INTERVALS_LONG_NAME = "cnv-intervals-vcf";
    public static final String CONTIG_PLOIDY_CALLS_LONG_NAME = "ploidy-calls-file";
    public static final String MIN_SIZE_LONG_NAME = "min-size";
    public static final String MAX_SIZE_LONG_NAME = "max-size";
    public static final String COPY_NEUTRAL_PRIOR_LONG_NAME = "copy-neutral-prior";

    @Argument(
            doc = "Germline copy number intervals VCF",
            fullName = COPY_NUMBER_INTERVALS_LONG_NAME
    )
    private File copyNumberPosteriorsFile;

    @Argument(
            doc = "Contig ploidy calls file. Can be specified for multiple samples.",
            fullName = CONTIG_PLOIDY_CALLS_LONG_NAME
    )
    private List<File> contigPloidyCallFiles;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFile;

    @Argument(
            doc = "Min event size",
            fullName = MIN_SIZE_LONG_NAME,
            minValue = 0,
            maxValue = Integer.MAX_VALUE,
            optional = true
    )
    private int minEventSize = 0;

    @Argument(
            doc = "Max event size",
            fullName = MAX_SIZE_LONG_NAME,
            minValue = 0,
            maxValue = Integer.MAX_VALUE,
            optional = true
    )
    private int maxEventSize = Integer.MAX_VALUE;

    @Argument(
            doc = "Prior probability of neutral copy number for " + COPY_NEUTRAL_PRIOR_BASIS_LENGTH +
                    "bp in regions where copy number calls are not available.",
            fullName = COPY_NEUTRAL_PRIOR_LONG_NAME,
            minValue = 0,
            maxValue = 1
    )
    private double copyNeutralPrior = 0.99;

    public final int COPY_NEUTRAL_PRIOR_BASIS_LENGTH = 1000;
    final String COPY_NUMBER_LOG_POSTERIORS_KEY = "CNLP";
    final String NEUTRAL_COPY_NUMBER_KEY = "NCN";

    private int numCopyStates;
    private List<IntegerCopyNumberState> copyStates;
    private VCFFileReader posteriorsReader;
    private GermlineCNVIntervalVariantDecoder cnvDecoder;
    private Map<String,Map<String,Integer>> sampleContigPloidyMap;
    private List<String> samples;
    private VariantContextWriter outputWriter;
    private String currentContig;
    private IntervalTree<Map<String,double[]>> currentTree;
    private SAMSequenceDictionary dictionary;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        posteriorsReader = new VCFFileReader(copyNumberPosteriorsFile);
        samples = getHeaderForVariants().getSampleNamesInOrder();
        final Collection<CalledContigPloidyCollection> contigPloidyCollections = contigPloidyCallFiles.stream()
                .map(SVCopyNumberPosteriors::readPloidyCalls)
                .collect(Collectors.toList());
        final Set<String> contigPloidySamples = contigPloidyCollections.stream()
                .map(p -> p.getMetadata().getSampleName())
                .collect(Collectors.toSet());
        final Set<String> svSamples = new HashSet<>(samples);
        final Set<String> cnvSamples = new HashSet<>(posteriorsReader.getFileHeader().getSampleNamesInOrder());
        if (!contigPloidySamples.containsAll(svSamples)) {
            final Set<String> missingSamples = new HashSet<>(svSamples);
            missingSamples.removeAll(contigPloidySamples);
            throw new UserException.BadInput("The following samples from the structural variant VCF was not in the contig ploidy calls: " + String.join(", ", missingSamples));
        }
        if (!cnvSamples.containsAll(svSamples)) {
            final Set<String> missingSamples = new HashSet<>(svSamples);
            missingSamples.removeAll(cnvSamples);
            throw new UserException.BadInput("The following samples from the structural variant VCF were not in the CNV VCF:" + String.join(", ", missingSamples));
        }
        cnvDecoder = new GermlineCNVIntervalVariantDecoder(contigPloidyCollections);
        sampleContigPloidyMap = cnvDecoder.getSampleContigPloidyMap();
        outputWriter = createVCFWriter(outputFile);
        outputWriter.writeHeader(composeHeader());
        currentContig = null;
        currentTree = null;
        initializeDefaultCopyStates();
    }

    private void initializeDefaultCopyStates() {
        final VariantContext exampleVariant = posteriorsReader.iterator().next();
        final SimpleInterval exampleInterval = new SimpleInterval(exampleVariant.getContig(), exampleVariant.getStart(), exampleVariant.getEnd());
        copyStates = cnvDecoder.parseGenotype(exampleVariant.getGenotype(0), exampleInterval)
                .getCopyNumberPosteriorDistribution().getIntegerCopyNumberStateList();
        numCopyStates = copyStates.size();
    }

    @Override
    public Object onTraversalSuccess() {
        outputWriter.close();
        return null;
    }

    private IntervalTree<Map<String,double[]>> getCurrentTree() {
        final SAMSequenceRecord contigRecord = dictionary.getSequence(currentContig);
        if (contigRecord == null) {
            throw new UserException.MissingContigInSequenceDictionary(currentContig, dictionary);
        }
        final Iterator<VariantContext> posteriorsIter = posteriorsReader.query(currentContig, 1, contigRecord.getSequenceLength());
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

    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (!isValidSize(variant, minEventSize, maxEventSize)) return;

        final String contig = variant.getContig();
        if (currentContig == null || contig != currentContig) {
            currentContig = contig;
            currentTree = getCurrentTree();
        }

        final SimpleInterval interval = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getEnd());
        final List<IntervalTree.Node<Map<String,double[]>>> posteriorsList = Lists.newArrayList(currentTree.overlappers(interval.getStart(), interval.getEnd()));
        final Map<String, CopyNumberPosteriorDistribution> variantPosteriors = samples.stream()
                .map(s -> new Tuple2<>(s, getEventPosterior(interval, posteriorsList, s)))
                .collect(Collectors.toMap(t -> t._1, t -> t._2));

        final VariantContextBuilder variantBuilder = new VariantContextBuilder(variant);
        final GenotypesContext genotypesContext = variantBuilder.getGenotypes();
        final List<Genotype> genotypesList = new ArrayList<>(samples.size());
        for (final Genotype genotype : genotypesContext) {
            final String sample = genotype.getSampleName();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);
            final CopyNumberPosteriorDistribution dist = variantPosteriors.get(sample);
            final List<Integer> copyStatePosterior = dist.getIntegerCopyNumberStateList().stream()
                    .map(s -> dist.getCopyNumberPosterior(s))
                    .map(p -> Integer.valueOf((int)(-10. * p / FastMath.log(10.))))
                    .collect(Collectors.toList());
            genotypeBuilder.attribute(COPY_NUMBER_LOG_POSTERIORS_KEY, copyStatePosterior);
            genotypeBuilder.attribute(NEUTRAL_COPY_NUMBER_KEY, getSamplePloidy(sample, contig));
            genotypesList.add(genotypeBuilder.make());
        }
        variantBuilder.genotypes(genotypesList);
        outputWriter.add(variantBuilder.make());
    }

    private VCFHeader composeHeader() {
        final VCFHeader variantHeader = getHeaderForVariants();
        final VCFHeader newHeader = new VCFHeader(variantHeader);
        getDefaultToolVCFHeaderLines().forEach(newHeader::addMetaDataLine);
        newHeader.addMetaDataLine(new VCFFormatHeaderLine(COPY_NUMBER_LOG_POSTERIORS_KEY, 1,
                VCFHeaderLineType.Integer, "Phred-scaled copy number posterior over the event region"));
        newHeader.addMetaDataLine(new VCFFormatHeaderLine(NEUTRAL_COPY_NUMBER_KEY, 1,
                VCFHeaderLineType.Integer, "Neutral copy number"));
        return newHeader;
    }

    private boolean isValidSize(final Locatable loc, final int min, final int max) {
        final int size = loc.getLengthOnReference();
        return size >= min && size <= max;
    }

    private int validatePosteriors(final LocatableCopyNumberPosteriorDistributionCollection posteriors) {
        if (posteriors.getRecords().isEmpty()) {
            throw new UserException.BadInput("Copy number posteriors collection for sample " + posteriors.getSampleName() + " is empty");
        }
        final int size = posteriors.getRecords().stream().map(LocatableCopyNumberPosteriorDistribution::getLengthOnReference)
                .collect(Collectors.groupingBy(i -> i)).entrySet().stream()
                .sorted(Comparator.comparingInt(e -> -e.getValue().size())).map(Map.Entry::getKey)
                .collect(Collectors.toList()).get(0);
        logger.info("Automatically determined interval size to be: " + size);
        final long nonUniformIntervals = posteriors.getRecords().stream().filter(r -> r.getLengthOnReference() != size).count();
        if (nonUniformIntervals > 0) {
            logger.warn("There were " + nonUniformIntervals + " copy number intervals for sample " + posteriors.getSampleName() + " not of size " + size);
        }
        return size;
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

    private static CalledContigPloidyCollection readPloidyCalls(final File file) {
        return new CalledContigPloidyCollection(file);
    }

    private CopyNumberPosteriorDistribution getEventPosterior(final SimpleInterval variantInterval,
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
}
