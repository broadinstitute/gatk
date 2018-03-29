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
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
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
public class DiscoverVariantsFromReadDepth extends GATKTool {

    private static final long serialVersionUID = 1L;

    public static final String BREAKPOINTS_LONG_NAME = "breakpoint-vcf";
    public static final String EVIDENCE_TARGET_LINKS_LONG_NAME = "evidence-target-links-file";
    public static final String SEGMENT_CALLS_LONG_NAME = "model-segments-file";
    public static final String COPY_RATIOS_LONG_NAME = "copy-ratio-file";
    public static final String ASSEMBLY_CONTIGS_LONG_NAME = "assembly-bam";

    @Argument(doc = "Breakpoints list (.vcf)", fullName = BREAKPOINTS_LONG_NAME)
    private String breakpointVCFPath;
    @Argument(doc = "Evidence target links file (.bedpe)", fullName = EVIDENCE_TARGET_LINKS_LONG_NAME)
    private String evidenceTargetLinksFilePath;
    @Argument(doc = "Model segments pipeline calls", fullName = SEGMENT_CALLS_LONG_NAME)
    private String segmentCallsFilePath;
    @Argument(doc = "Copy number ratio file (.tsv or .hdf5)", fullName = COPY_RATIOS_LONG_NAME)
    private String copyRatioFilePath;
    @Argument(doc = "Assembly contigs (.bam or .sam)", fullName = ASSEMBLY_CONTIGS_LONG_NAME, optional = true)
    private String assemblyBamPath;
    @Argument(doc = "Output intervals path (.bed)", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputPath;

    private final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepth arguments = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepth();
    private LargeSimpleSVCaller largeSimpleSVCaller;
    private SAMSequenceDictionary dictionary;
    private Collection<LargeSimpleSV> events;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        logger.info("Loading assembly...");
        final Collection<GATKRead> assembly = getReads(assemblyBamPath);
        logger.info("Loading breakpoints...");
        final Collection<VariantContext> breakpoints = readVCF(breakpointVCFPath, dictionary);
        logger.info("Loading copy ratio segment calls...");
        final CalledCopyRatioSegmentCollection copyRatioSegments = getCalledCopyRatioSegments(segmentCallsFilePath);
        logger.info("Loading evidence target links...");
        final Collection<EvidenceTargetLink> evidenceTargetLinks = getEvidenceTargetLinks(evidenceTargetLinksFilePath, dictionary);
        logger.info("Loading copy ratio bins...");
        final CopyRatioCollection copyRatios = getCopyRatios(copyRatioFilePath);

        largeSimpleSVCaller = new LargeSimpleSVCaller(breakpoints, assembly, evidenceTargetLinks, copyRatios, copyRatioSegments, dictionary, arguments);
    }

    @Override
    public void traverse() {
        events = largeSimpleSVCaller.callEvents(progressMeter);
    }

    @Override
    public Object onTraversalSuccess() {
        writeEvents(outputPath, events, dictionary);
        return null;
    }

    /**
     * Writes events collection to BED file
     */
    private void writeEvents(final String outputPath, final Collection<LargeSimpleSV> events, final SAMSequenceDictionary dictionary) {
        try (final OutputStream outputStream = BucketUtils.createFile(outputPath)) {
            outputStream.write((LargeSimpleSV.getBedHeader() + "\n").getBytes());
            for (final LargeSimpleSV event : events) {
                outputStream.write((event.toBedString(dictionary, arguments.counterEvidencePseudocount) + "\n").getBytes());
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

    private static CalledCopyRatioSegmentCollection getCalledCopyRatioSegments(final String path) {
        final File modelSegmentsCallsFile = new File(path);
        return new CalledCopyRatioSegmentCollection(modelSegmentsCallsFile);
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
