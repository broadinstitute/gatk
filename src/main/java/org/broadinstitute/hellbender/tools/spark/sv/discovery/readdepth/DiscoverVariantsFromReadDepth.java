package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.ModelSegments;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Call structural variants from the output of the SV and Model Segments CNV pipelines
 *
 * <p>See StructuralVariationDiscoveryPipelineSpark and ModelSegments for information on these pipelines.</p>
 *
 * <p>This tool currently only calls large (> 1 kbp) tandem duplications.</p>
 *
 * <h3>Methods</h3>
 *
 * <p>The tool scans the genome to look for events supported by a combination of local assembly contigs, breakpoint pairs,
 * target evidence links, copy number ratio, and copy ratio segment calls.</p>
 *
 * <p>To detect tandem duplications, putative event intervals are generated using intrachromosomal breakpoint (BND)
 * pairs and groups of discordant read pairs and split reads (target evidence links). All target evidence links overlapping
 * the current interval are retrieved and evaluated as to whether they are consistent with the presence of a tandem duplication over
 * the interval. In particular, the strands of the target evidence links must be -/+ (i.e. "outie" read pair orientation).
 * Target evidence links that are of different strands (+/-, -/-, +/+) or do not overlap the ends of the interval within
 * some tolerance are marked as "counter-evidence" that suggest either a different kind of event or an overlapping event
 * is occurring on the interval (a pseudocount is added if there is no evidence). Intervals where the ratio of evidence
 * to counter-evidence are below a threshold are filtered out.</p>
 *
 * <p>Intervals with sufficient evidence to counter-evidence ratio are then tested against calls from the Model Segments
 * CNV pipeline. Intervals are called if they have at least 50% reciprocal overlap with an amplification call.</p>
 *
 * <p>If there is no matching Model Segments call, the Viterbi algorithm is run on the provided copy number ratio bins
 * using an HMM of copy number modeled with Gaussian distributions. The interval is then called if the most likely
 * sequence of states consists of a sufficient number of elevated copy number states (> 2).</p>
 *
 * <h3>Input</h3>
 *
 * <ul>
 *     <li>Breakpoint VCF from StructuralVariationDiscoveryPipelineSpark</li>
 *     <li>Evidence target links .bedpe file from StructuralVariationDiscoveryPipelineSpark</li>
 *     <li>ModelSegments .seg calls file</li>
 *     <li>Copy number ratio .tsv file</li>
 *     <li>SAM/BAM of aligned assembly contigs from StructuralVariationDiscoveryPipelineSpark</li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <p>Bed file containing the following columns</p>
 *
 * <ul>
 *      <li>contig name</li>
 *      <li>start position</li>
 *      <li>end position</li>
 *      <li>Read pair evidence</li>
 *      <li>Split read evidence</li>
 *      <li>Read pair counter-evidence</li>
 *      <li>Split read counter-evidence</li>
 *      <li>Ratio score</li>
 *      <li>Associated breakpoint pair (if any)</li>
 * </ul>
 *
 * <p>The tool may also optionally produce:</p>
 *
 * <ul>
 *     <li>BAM file of all reads annotated with the NCBI taxonomy IDs of mapped organisms</li>
 *     <li>Metrics file with the number of mapped and unmapped reads</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 * gatk DiscoverVariantsFromReadDepth  \
 *   --breakpoint-vcf breakpoints.vcf \
 *   --evidence-target-links-file target_links.bedpe \
 *   --model-segments-file segments.seg \
 *   --copy-ratio-file copy_ratios.tsv \
 *   --assembly-bam assemblies.sam \
 *   --output output.bed
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@ExperimentalFeature
@CommandLineProgramProperties(
        oneLineSummary = "Finds large tandem duplications",
        summary = "This tool uses the output from StructuralVariationDiscoveryPipelineSpark and the Model Segments CNV pipeline " +
                "to call large tandem duplications.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class DiscoverVariantsFromReadDepth extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    public static final String SV_CALLS_LONG_NAME = "sv-calls-vcf";
    public static final String BREAKPOINTS_LONG_NAME = "breakpoint-vcf";
    public static final String EVIDENCE_TARGET_LINKS_LONG_NAME = "evidence-target-links-file";
    public static final String CNV_CALLS_LONG_NAME = "cnv-segments-file";
    public static final String COPY_RATIOS_LONG_NAME = "copy-ratio-file";
    public static final String ASSEMBLY_CONTIGS_LONG_NAME = "assembly-bam";
    public static final String FILTERED_CALLS_LONG_NAME = "filtered-calls";
    public static final String HIGH_COVERAGE_INTERVALS_LONG_NAME = "high-coverage-intervals";

    @Argument(doc = "Breakpoints list (.vcf)", fullName = BREAKPOINTS_LONG_NAME)
    private String breakpointVCFPath;
    @Argument(doc = "Structural variant calls list (.vcf)", fullName = SV_CALLS_LONG_NAME)
    private String svCallVCFPath;
    @Argument(doc = "Evidence target links file (.bedpe)", fullName = EVIDENCE_TARGET_LINKS_LONG_NAME)
    private String evidenceTargetLinksFilePath;
    @Argument(doc = "CNV segments calls (gCNV .vcf or ModelSegments .seg)", fullName = CNV_CALLS_LONG_NAME)
    private String segmentCallsFilePath;
    @Argument(doc = "Copy number ratio file (.tsv or .hdf5)", fullName = COPY_RATIOS_LONG_NAME)
    private String copyRatioFilePath;
    @Argument(doc = "Assembly contigs (.bam or .sam)", fullName = ASSEMBLY_CONTIGS_LONG_NAME, optional = true)
    private String assemblyBamPath;
    @Argument(doc = "Output intervals path (.bed)", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputPath;
    @Argument(doc = "Output filtered calls path (.vcf)",fullName = FILTERED_CALLS_LONG_NAME)
    private String filteredCallsPath;
    @Argument(doc = "High coverage intervals file path",fullName = HIGH_COVERAGE_INTERVALS_LONG_NAME)
    private String highCoverageIntervalsPath;
    @Argument(doc = "Mappable intervals file path",fullName = "mappable-intervals")
    private String mappableIntervalsPath;
    @Argument(doc = "One or more blacklisted intervals files",fullName = "blacklist")
    private List<String > blacklistedPaths;
    @Argument(doc = "Sequence dictionary",fullName = "sequence-dictionary")
    private String sequenceDictionaryPath;
    @Argument(doc = "Truth set (.vcf)",fullName = "truth-intervals", optional = true)
    private String truthIntervalsPath;

    private final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection();
    private LargeSimpleSVCaller largeSimpleSVCaller;
    private SAMSequenceDictionary dictionary;
    private Collection<ReadDepthEvent> events;
    private List<VariantContext> filteredCalls;

    @Override
    public void runTool(final JavaSparkContext ctx) {
        //ctx.hadoopConfiguration().set("fs.gs.project.id", "broad-dsde-methods");
        //ctx.hadoopConfiguration().set("google.cloud.auth.service.account.json.keyfile", "/Users/markw/IdeaProjects/gatk/DSDE-Methods-2fcf61fced3a.json");
        dictionary = ReferenceUtils.loadFastaDictionary(BucketUtils.openFile(sequenceDictionaryPath));
        final JavaRDD<GATKRead> reads = getUnfilteredReads();
        logger.info("Loading assembly...");
        final Collection<GATKRead> assembly = getReads(assemblyBamPath);
        logger.info("Loading high coverage intervals...");
        final List<SVInterval> highCoverageIntervals = getIntervals(highCoverageIntervalsPath, dictionary);
        logger.info("Loading high coverage intervals...");
        final List<SVInterval> mappableIntervals = getIntervals(mappableIntervalsPath, dictionary);
        final List<SVInterval> blacklist = new ArrayList<>();
        for (final String path : blacklistedPaths) {
            logger.info("Loading blacklist at " + path);
            blacklist.addAll(getIntervals(path, dictionary));
        }
        final Collection<VariantContext> truthSet;
        if (truthIntervalsPath == null) {
            truthSet = null;
        } else {
            logger.info("Loading truth set...");
            truthSet = readVCF(truthIntervalsPath, dictionary);
        }
        logger.info("Loading breakpoints...");
        final Collection<VariantContext> breakpointCalls = readVCF(breakpointVCFPath, dictionary);
        logger.info("Loading SV calls...");
        final Collection<VariantContext> svCalls = readVCF(svCallVCFPath, dictionary);
        logger.info("Loading copy ratio segment calls...");
        final CalledCopyRatioSegmentCollection copyRatioSegments = getCalledCopyRatioSegments(segmentCallsFilePath, dictionary);
        logger.info("Loading evidence target links...");
        final Collection<EvidenceTargetLink> evidenceTargetLinks = getEvidenceTargetLinks(evidenceTargetLinksFilePath, dictionary);
        logger.info("Loading copy ratio bins...");
        final CopyRatioCollection copyRatios = getCopyRatios(copyRatioFilePath);

        final List<LargeSimpleSV> largeSimpleSVs = Collections.emptyList(); //TODO

        largeSimpleSVCaller = new LargeSimpleSVCaller(reads, largeSimpleSVs, svCalls, truthSet, assembly, evidenceTargetLinks, copyRatios, copyRatioSegments, highCoverageIntervals, mappableIntervals, blacklist, dictionary, arguments);

        final Tuple2<Collection<ReadDepthEvent>,List<VariantContext>> result = largeSimpleSVCaller.callEvents(readArguments.getReadFilesNames().get(0), ctx, null);
        events = result._1;
        filteredCalls = result._2;

        writeEvents(outputPath, events, dictionary);
        SVVCFWriter.writeVCF(filteredCalls, filteredCallsPath, dictionary, logger);
    }

    /**
     * Writes events collection to BED file
     */
    private void writeEvents(final String outputPath, final Collection<ReadDepthEvent> events, final SAMSequenceDictionary dictionary) {
        try (final OutputStream outputStream = BucketUtils.createFile(outputPath)) {
            outputStream.write(("#" + ReadDepthEvent.getBedHeader() + "\n").getBytes());
            for (final ReadDepthEvent event : events) {
                outputStream.write((event.toBedString(dictionary) + "\n").getBytes());
            }
        } catch (final IOException e) {
            throw new GATKException("Error writing output BED file", e);
        }
    }

    private static Collection<VariantContext> readVCF(final String vcfPath, final SAMSequenceDictionary dictionary) {
        try (final FeatureDataSource<VariantContext> dataSource =
                     new FeatureDataSource<>(vcfPath, null, 0, VariantContext.class)) {
            Collection<VariantContext> variants = Utils.stream(dataSource.iterator()).collect(Collectors.toList());
            if (variants.stream().anyMatch(variant -> dictionary.getSequenceIndex(variant.getContig()) < 0)) {
                throw new UserException("Found a VCF entry contig that does not appear in the dictionary.");
            }
            return variants;
        }
    }

    private static List<SVInterval> getIntervals(final String path, final SAMSequenceDictionary dictionary) {
        if (path != null) {
            final GenomeLocParser genomeLocParser = new GenomeLocParser(dictionary);
            return IntervalUtils.parseIntervalArguments(genomeLocParser, path).stream()
                    .map(genomeLoc -> new SVInterval(genomeLoc.getContigIndex(), genomeLoc.getStart(), genomeLoc.getEnd()))
                    .collect(Collectors.toList());
        }
        return Collections.emptyList();
    }

    private static Collection<GATKRead> getReads(final String bamPath) {
        if (bamPath != null) {
            final ReadsDataSource reads = new ReadsDataSource(IOUtils.getPath(bamPath));
            return Utils.stream(reads.iterator()).collect(Collectors.toList());
        }
        return Collections.emptyList();
    }

    private static CopyRatioCollection getCopyRatios(final String path) {
        final File copyRatioFile = new File(path);
        return new CopyRatioCollection(copyRatioFile);
    }

    private static CalledCopyRatioSegmentCollection getCalledCopyRatioSegments(final String path, final SAMSequenceDictionary dictionary) {
        final File file = new File(path);
        if (path.toLowerCase().endsWith(ModelSegments.SEGMENTS_FILE_SUFFIX.toLowerCase())) {
            return new CalledCopyRatioSegmentCollection(file);
        } else if (path.toLowerCase().endsWith(IOUtil.VCF_FILE_EXTENSION.toLowerCase())) {
            final VCFFileReader reader = new VCFFileReader(file, false);
            final List<CalledCopyRatioSegment> segments = Utils.stream(reader.iterator()).map(variantContext ->  {
                final SimpleInterval interval = new SimpleInterval(variantContext.getContig(), variantContext.getStart(), variantContext.getEnd());
                final GenotypesContext genotypesContext = variantContext.getGenotypes();
                if (genotypesContext.isEmpty()) {
                    throw new UserException.BadInput("No genotypes found in variant context " + variantContext.getID());
                }
                final Genotype genotype = genotypesContext.get(0);
                if (!genotype.hasExtendedAttribute("NP")) {
                    throw new UserException.BadInput("Number of points genotype not found in variant context " + variantContext.getID());
                }
                if (!genotype.hasExtendedAttribute("CN")) {
                    throw new UserException.BadInput("Copy number genotype not found in variant context " + variantContext.getID());
                }
                final int copyNumber = Integer.valueOf((String) genotype.getExtendedAttribute("CN"));
                final int numPoints = Integer.valueOf((String) genotype.getExtendedAttribute("NP"));
                final CopyRatioSegment segment = new CopyRatioSegment(interval, numPoints, Math.log(copyNumber) / Math.log(2.0));
                final CalledCopyRatioSegment.Call call;
                if (copyNumber == 2) {
                    call = CalledCopyRatioSegment.Call.NEUTRAL;
                } else if (copyNumber > 2) {
                    call = CalledCopyRatioSegment.Call.AMPLIFICATION;
                } else {
                    call = CalledCopyRatioSegment.Call.DELETION;
                }
                return new CalledCopyRatioSegment(segment, call);
            }).collect(Collectors.toList());
            final String sampleName = reader.getFileHeader().getSampleNamesInOrder().get(0);
            return new CalledCopyRatioSegmentCollection(new SimpleSampleLocatableMetadata(sampleName, dictionary), segments);
        }
        throw new UserException.BadInput("CNV segments file must be " + ModelSegments.SEGMENTS_FILE_SUFFIX + " or " + IOUtil.VCF_FILE_EXTENSION);
    }

    private static Collection<EvidenceTargetLink> getEvidenceTargetLinks(final String path, final SAMSequenceDictionary dictionary) {
        final File file = new File(path);
        final Collection<EvidenceTargetLink> links = new ArrayList<>();
        final String[] evidenceTargetFileLines = new String(IOUtils.readFileIntoByteArray(file)).split("\n");
        for (final String line : evidenceTargetFileLines) {
            links.add(EvidenceTargetLink.fromBedpeString(line.trim(), dictionary));
        }
        return links;
    }

}
