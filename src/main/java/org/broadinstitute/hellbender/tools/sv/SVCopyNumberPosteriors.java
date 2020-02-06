package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
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
    public static final String COPY_NEUTRAL_PRIOR_LONG_NAME = "copy-neutral-prior";
    public static final String CNV_BND_PROB_THRESHOLD_LONG_NAME = "cnv-bnd-prob-thresh";
    public static final String CNV_BND_SAMPLE_THRESHOLD_LONG_NAME = "cnv-bnd-sample-thresh";

    @Argument(
            doc = "Germline copy number intervals VCF. Can be specified more than once for runs with different bin sizes.",
            fullName = COPY_NUMBER_INTERVALS_LONG_NAME
    )
    private List<File> copyNumberPosteriorsFiles;

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
            doc = "Prior probability of neutral copy number for " + COPY_NEUTRAL_PRIOR_BASIS_LENGTH +
                    "bp in regions where copy number calls are not available.",
            fullName = COPY_NEUTRAL_PRIOR_LONG_NAME,
            minValue = 0,
            maxValue = 1
    )
    private double copyNeutralPrior = 0.99;

    @Argument(
            doc = "Threshold phred-scaled alt copy state probability for BND conversion of CNVs.",
            fullName = CNV_BND_PROB_THRESHOLD_LONG_NAME,
            minValue = 0
    )
    private int cnvBndProbThresh = 20;

    @Argument(
            doc = "Threshold carrier sample fraction meeting probability threshold for BND conversion of CNVs.",
            fullName = CNV_BND_SAMPLE_THRESHOLD_LONG_NAME,
            minValue = 0.,
            maxValue = 1.
    )
    private double cnvBndSampleThresh = 0.5;

    public final int COPY_NEUTRAL_PRIOR_BASIS_LENGTH = 1000;
    final String COPY_NUMBER_LOG_POSTERIORS_KEY = "CNLP";
    final String NEUTRAL_COPY_NUMBER_KEY = "NCN";
    public final List<StructuralVariantType> CNV_TYPES = Lists.newArrayList(StructuralVariantType.DEL, StructuralVariantType.DUP);

    private final Map<String,IntervalTree<Object>> whitelistedIntervalTreeMap = new HashMap<>();
    private int numCopyStates;
    private List<IntegerCopyNumberState> copyStates;
    private List<VCFFileReader> posteriorsReaders;
    private GermlineCNVIntervalVariantDecoder cnvDecoder;
    private Map<String,Map<String,Integer>> sampleContigPloidyMap;
    private List<String> samples;
    private VariantContextWriter outputWriter;
    private String currentContig;
    private List<IntervalTree<Map<String,double[]>>> currentPosteriorsTreeList;
    private SAMSequenceDictionary dictionary;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        posteriorsReaders = copyNumberPosteriorsFiles.stream().map(VCFFileReader::new).collect(Collectors.toList());
        samples = getHeaderForVariants().getSampleNamesInOrder();
        final Collection<CalledContigPloidyCollection> contigPloidyCollections = contigPloidyCallFiles.stream()
                .map(SVCopyNumberPosteriors::readPloidyCalls)
                .collect(Collectors.toList());
        final Set<String> contigPloidySamples = contigPloidyCollections.stream()
                .map(p -> p.getMetadata().getSampleName())
                .collect(Collectors.toSet());
        final Set<String> svSamples = new HashSet<>(samples);
        final Set<String> cnvSamples = getCNVSamples();
        if (!contigPloidySamples.containsAll(svSamples)) {
            final Set<String> missingSamples = new HashSet<>(svSamples);
            missingSamples.removeAll(contigPloidySamples);
            throw new UserException.BadInput("The following samples from the structural variant VCF were not in the contig ploidy calls: " + String.join(", ", missingSamples));
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
        currentPosteriorsTreeList = new ArrayList<>(copyNumberPosteriorsFiles.size());
        initializeDefaultCopyStates();
        loadIntervalTree();
    }

    private Set<String> getCNVSamples() {
        final Set<String> cnvSamplesSet = new HashSet<>(posteriorsReaders.get(0).getFileHeader().getSampleNamesInOrder());
        for (int i = 1; i < posteriorsReaders.size(); i++) {
            final Set<String> set = new HashSet<>(posteriorsReaders.get(i).getFileHeader().getSampleNamesInOrder());
            if (!set.equals(cnvSamplesSet)) {
                throw new UserException.BadInput("CNV VCFs do not contain identical samples.");
            }
        }
        return Collections.unmodifiableSet(cnvSamplesSet);
    }

    private void loadIntervalTree() {
        if (getTraversalIntervals() == null) {
            for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
                whitelistedIntervalTreeMap.put(sequence.getSequenceName(), new IntervalTree<>());
                whitelistedIntervalTreeMap.get(sequence.getSequenceName()).put(1, sequence.getSequenceLength(), null);
            }
        } else {
            for (final SimpleInterval interval : getTraversalIntervals()) {
                whitelistedIntervalTreeMap.putIfAbsent(interval.getContig(), new IntervalTree<>());
                whitelistedIntervalTreeMap.get(interval.getContig()).put(interval.getStart(), interval.getEnd(), null);
            }
        }
    }

    private void initializeDefaultCopyStates() {
        final VariantContext exampleVariant = posteriorsReaders.get(0).iterator().next();
        copyStates = getCopyNumberStates(exampleVariant);
        numCopyStates = copyStates.size();
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

    @Override
    public Object onTraversalSuccess() {
        outputWriter.close();
        return null;
    }

    private List<IntervalTree<Map<String,double[]>>> getCurrentPosteriorsTreeList() {
        return posteriorsReaders.stream().map(this::getCurrentPosteriorsTree).collect(Collectors.toList());
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

    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final SVCallRecord call = SVCallRecordWithEvidence.create(variant);
        if (!SVCluster.isValidSize(call, minEventSize) || !SVCluster.isWhitelisted(call, whitelistedIntervalTreeMap)) return;

        final String contig = variant.getContig();
        if (contig != currentContig) {
            currentContig = contig;
            currentPosteriorsTreeList = getCurrentPosteriorsTreeList();
        }

        final Map<String, CopyNumberPosteriorDistribution> variantPosteriors = samples.stream()
                .map(s -> new Tuple2<>(s, getCallPosterior(call, s)))
                .collect(Collectors.toMap(t -> t._1, t -> t._2));
        outputWriter.add(makeVariant(variant, variantPosteriors));
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

    private CopyNumberPosteriorDistribution getCallPosterior(final SVCallRecord call,
                                                             final String sample) {
        if (!CNV_TYPES.contains(call.getType())) {
            return getDefaultPosterior();
        }
        if (!call.getContig().equals(call.getEndContig())) {
            throw new IllegalArgumentException("Encountered interchromosomal CNV call");
        }
        final SimpleInterval interval = new SimpleInterval(call.getContig(), call.getStart(), call.getEnd());
        final List<IntervalTree.Node<Map<String,double[]>>> posteriorsList = getBestPosteriors(interval);
        return getIntervalPosterior(interval, posteriorsList, sample);
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

    private CopyNumberPosteriorDistribution getDefaultPosterior() {
        final double p = FastMath.log(1. / numCopyStates);
        final Map<IntegerCopyNumberState,Double> dist = new HashMap<>(SVUtils.hashMapCapacity(numCopyStates));
        for (final IntegerCopyNumberState s : copyStates) {
            dist.put(s, p);
        }
        return new CopyNumberPosteriorDistribution(dist);
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

    private VariantContext makeVariant(final VariantContext variant,
                                       final Map<String, CopyNumberPosteriorDistribution> variantPosteriors) {
        final VariantContextBuilder variantBuilder = new VariantContextBuilder(variant);
        final List<Genotype> genotypesList = new ArrayList<>(samples.size());
        final StructuralVariantType svType = variant.getStructuralVariantType();
        final Map<String, List<Integer>> phredScaledPosteriors = getPhredScaledPosteriors(variantPosteriors);

        // Convert CNVs lacking depth support to BNDs
        final Allele altAllele = Allele.create("<" + StructuralVariantType.BND + ">", false);
        final Allele refAllele = Allele.REF_N;
        final List<String> algorithms = variant.getAttributeAsStringList(SVCluster.ALGORITHMS_ATTRIBUTE, null);
        final boolean convertToBnd = (svType.equals(StructuralVariantType.DEL) || svType.equals(StructuralVariantType.DUP))
                && !(algorithms.size() == 1 && algorithms.get(0).equals(SVCluster.DEPTH_ALGORITHM))
                && !hasDepthSupport(variant, phredScaledPosteriors);
        if (convertToBnd) {
            variantBuilder.alleles(Lists.newArrayList(refAllele, altAllele));
            variantBuilder.attribute(VCFConstants.SVTYPE, StructuralVariantType.BND);
            variantBuilder.attribute(SVCluster.SVLEN_ATTRIBUTE, -1);
        }

        final GenotypesContext genotypesContext = variantBuilder.getGenotypes();
        for (final Genotype genotype : genotypesContext) {
            final String sample = genotype.getSampleName();
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(genotype);
            if (convertToBnd && !genotype.isHomRef()) {
                genotypeBuilder.alleles(Lists.newArrayList(refAllele, altAllele));
            }
            genotypeBuilder.attribute(COPY_NUMBER_LOG_POSTERIORS_KEY, phredScaledPosteriors.get(sample));
            genotypeBuilder.attribute(NEUTRAL_COPY_NUMBER_KEY, getSamplePloidy(sample, variant.getContig()));
            genotypesList.add(genotypeBuilder.make());
        }
        variantBuilder.genotypes(genotypesList);
        return variantBuilder.make();
    }

    private Map<String, List<Integer>> getPhredScaledPosteriors(final Map<String, CopyNumberPosteriorDistribution> variantPosteriors) {
        final Map<String, List<Integer>> phredScaledPosteriors = new HashMap<>(SVUtils.hashMapCapacity(variantPosteriors.size()));
        for (final String sample : variantPosteriors.keySet()) {
            final CopyNumberPosteriorDistribution dist = variantPosteriors.get(sample);
            final List<Integer> copyStatePosterior = dist.getIntegerCopyNumberStateList().stream()
                    .map(s -> dist.getCopyNumberPosterior(s))
                    .map(p -> Integer.valueOf((int)(-10. * p / FastMath.log(10.))))
                    .collect(Collectors.toList());
            phredScaledPosteriors.put(sample, copyStatePosterior);
        }
        return phredScaledPosteriors;
    }

    private boolean hasDepthSupport(final VariantContext variant,
                                    final Map<String, List<Integer>> phredScaledPosteriors) {
        final StructuralVariantType svType = variant.getStructuralVariantType();
        if (!svType.equals(StructuralVariantType.DEL) && !svType.equals(StructuralVariantType.DUP)) {
            throw new IllegalArgumentException("Variant was not a CNV");
        }
        int samplesWithDepthSupport = 0;
        for (final String sample : phredScaledPosteriors.keySet()) {
            final List<Integer> copyStatePosterior = phredScaledPosteriors.get(sample);
            final int samplePloidy = getSamplePloidy(sample, variant.getContig());
            if (svType.equals(StructuralVariantType.DEL)) {
                for (int i = 0; i < samplePloidy; i++) {
                    if (copyStatePosterior.get(i) <= cnvBndProbThresh) {
                        samplesWithDepthSupport++;
                        break;
                    }
                }
            } else {
                for (int i = samplePloidy + 1; i < copyStatePosterior.size(); i++) {
                    if (copyStatePosterior.get(i) <= cnvBndProbThresh) {
                        samplesWithDepthSupport++;
                        break;
                    }
                }
            }
        }
        final int numCarriers = variant.getNSamples() - (int) variant.getGenotypes().stream().filter(g -> g.isHomRef()).count();
        final int minSamples = (int) Math.ceil(numCarriers * cnvBndSampleThresh);
        return samplesWithDepthSupport >= minSamples;
    }
}
