package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Paths;
import java.util.*;
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
public class StructuralVariantIntervalsForCNV extends GATKTool {

    private static final long serialVersionUID = 1L;

    public static final String SV_CALLS_LONG_NAME = "sv-calls-vcf";
    public static final String BREAKPOINTS_LONG_NAME = "breakpoint-vcf";
    public static final String EVIDENCE_TARGET_LINKS_LONG_NAME = "evidence-target-links-file";
    public static final String ASSEMBLY_CONTIGS_LONG_NAME = "assembly-bam";
    public static final String HIGH_COVERAGE_INTERVALS_LONG_NAME = "high-coverage-intervals";

    @Argument(doc = "Breakpoints list (.vcf)", fullName = BREAKPOINTS_LONG_NAME)
    private String breakpointVCFPath;
    @Argument(doc = "Structural variant calls list (.vcf)", fullName = SV_CALLS_LONG_NAME)
    private String svCallVCFPath;
    @Argument(doc = "Evidence target links file (.bedpe)", fullName = EVIDENCE_TARGET_LINKS_LONG_NAME)
    private String evidenceTargetLinksFilePath;
    @Argument(doc = "Assembly contigs (.bam or .sam)", fullName = ASSEMBLY_CONTIGS_LONG_NAME, optional = true)
    private String assemblyBamPath;
    @Argument(doc = "Output directory", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputPath;
    @Argument(doc = "Sample name for output",fullName = "sample-name")
    private String sampleName;
    @Argument(doc = "High coverage intervals file path",fullName = HIGH_COVERAGE_INTERVALS_LONG_NAME)
    private String highCoverageIntervalsPath;
    @Argument(doc = "Mappable intervals file path",fullName = "mappable-intervals")
    private String mappableIntervalsPath;
    @Argument(doc = "One or more blacklisted intervals files",fullName = "blacklist")
    private List<String> blacklistedPaths;

    private final StructuralVariationDiscoveryArgumentCollection.StructuralVariantIntervalsForCNV arguments = new StructuralVariationDiscoveryArgumentCollection.StructuralVariantIntervalsForCNV();
    private StructuralVariantIntervalFinder structuralVariantIntervalFinder;
    private SAMSequenceDictionary dictionary;

    @Override
    public void traverse() {
        //ctx.hadoopConfiguration().set("fs.gs.project.id", "broad-dsde-methods");
        //ctx.hadoopConfiguration().set("google.cloud.auth.service.account.json.keyfile", "/Users/markw/IdeaProjects/gatk/DSDE-Methods-2fcf61fced3a.json");
        dictionary = getBestAvailableSequenceDictionary();
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
        logger.info("Loading breakpoints...");
        final Collection<VariantContext> breakpointCalls = readVCF(breakpointVCFPath, dictionary);
        logger.info("Loading SV calls...");
        final Collection<VariantContext> svCalls = readVCF(svCallVCFPath, dictionary);
        logger.info("Loading evidence target links...");
        final Collection<EvidenceTargetLink> evidenceTargetLinks = getEvidenceTargetLinks(evidenceTargetLinksFilePath, dictionary);

        structuralVariantIntervalFinder = new StructuralVariantIntervalFinder(breakpointCalls, svCalls, assembly, evidenceTargetLinks, highCoverageIntervals, mappableIntervals, blacklist, dictionary, arguments);
        final Map<String,List<TypedSVInterval>> svIntervals = structuralVariantIntervalFinder.getIntervals(progressMeter);
        writeEvents(outputPath, svIntervals);
    }

    /**
     * Writes events collection to BED file
     */
    private void writeEvents(final String outputDirectory, final Map<String,List<TypedSVInterval>> svIntervals) {
        for (final String key : svIntervals.keySet()) {
            final String fileName = sampleName + "-" + key + ".bed";
            final String outputPath = Paths.get(outputDirectory, fileName).toAbsolutePath().toString();
            try (final OutputStream outputStream = BucketUtils.createFile(outputPath)) {
                outputStream.write(("#CHR\tPOS\tEND\tTYPE\n").getBytes());
                for (final TypedSVInterval interval : svIntervals.get(key)) {
                    outputStream.write((dictionary.getSequence(interval.getContig()).getSequenceName() + "\t" + interval.getStart() + "\t" + interval.getEnd() + "\t" + interval.getType().name() + "\n").getBytes());
                }
            } catch (final IOException e) {
                throw new GATKException("Error writing output to " + outputPath, e);
            }
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
