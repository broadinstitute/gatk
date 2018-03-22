package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
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
 * gatk StructuralVariationDiscoveryWithDepth  \
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
public class StructuralVariationDiscoveryWithDepth extends GATKTool {

    private static final long serialVersionUID = 1L;
    private static final int HIGH_DEPTH_COVERAGE_PEAK_FACTOR = 7;

    @Argument(doc = "Breakpoint VCF", fullName = "breakpoint-vcf")
    private String breakpointVCFPath;
    @Argument(doc = "Evidence target links bedpe file", fullName = "evidence-target-links-file")
    private String evidenceTargetLinksFilePath;
    @Argument(doc = "Model segments pipeline calls", fullName = "model-segments-file")
    private String segmentCallsFilePath;
    @Argument(doc = "Copy number ratio tsv file", fullName = "copy-ratio-file")
    private String copyRatioFilePath;
    @Argument(doc = "Assembly SAM/BAM", fullName = "assembly-bam", optional = true)
    private String assemblyBamPath;
    @Argument(doc = "Output path", fullName = "output")
    private String outputPath;

    final StructuralVariationDiscoveryArgumentCollection.SimpleVariantDiscoveryArgumentCollection arguments = new StructuralVariationDiscoveryArgumentCollection.SimpleVariantDiscoveryArgumentCollection(HIGH_DEPTH_COVERAGE_PEAK_FACTOR);

    public void traverse() {
        final SAMSequenceDictionary dictionary = getBestAvailableSequenceDictionary();
        logger.info("Loading assembly...");
        final List<GATKRead> assembly = getReads(assemblyBamPath);
        logger.info("Loading breakpoints...");
        final List<VariantContext> breakpoints = readVCF(breakpointVCFPath, dictionary);
        logger.info("Loading copy ratio segment calls...");
        final CalledCopyRatioSegmentCollection copyRatioSegments = getCalledCopyRatioSegments(segmentCallsFilePath);
        logger.info("Loading evidence target links...");
        final List<EvidenceTargetLink> evidenceTargetLinks = getEvidenceTargetLinks(evidenceTargetLinksFilePath, dictionary);
        logger.info("Loading copy ratio bins...");
        final CopyRatioCollection copyRatios = getCopyRatios(copyRatioFilePath);

        logger.info("Working...");
        final DepthBasedSVCaller depthBasedSVCaller = new DepthBasedSVCaller(breakpoints, assembly, evidenceTargetLinks, copyRatios, copyRatioSegments, dictionary, arguments);
        final List<LargeSimpleSV> events = depthBasedSVCaller.callEvents();
        writeEvents(outputPath, events, dictionary);
    }

    private void writeEvents(final String outputPath, final Collection<LargeSimpleSV> events, final SAMSequenceDictionary dictionary) {
        try (final OutputStream outputStream = BucketUtils.createFile(outputPath)) {
            for (final LargeSimpleSV event : events) {
                outputStream.write((event.getString(dictionary, arguments.COUNTEREVIDENCE_PSEUDOCOUNT) + "\n").getBytes());
            }
        } catch (final IOException e) {
            throw new GATKException("Error writing output BED file", e);
        }
    }

    private static List<VariantContext> readVCF(final String vcfPath, final SAMSequenceDictionary dictionary) {
        try (final FeatureDataSource<VariantContext> dataSource =
                     new FeatureDataSource<>(vcfPath, null, 0, VariantContext.class)) {
            List<VariantContext> variants = Utils.stream(dataSource.iterator()).collect(Collectors.toList());
            if (variants.stream().anyMatch(variant -> dictionary.getSequenceIndex(variant.getContig()) < 0)) {
                throw new UserException("Found a VCF entry contig that does not appear in the dictionary.");
            }
            return variants;
        }
    }

    private static List<GATKRead> getReads(final String bamPath) {
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

    private static List<EvidenceTargetLink> getEvidenceTargetLinks(final String path, final SAMSequenceDictionary dictionary) {
        final File file = new File(path);
        final List<EvidenceTargetLink> links = new ArrayList<>();
        final String[] evidenceTargetFileLines = new String(IOUtils.readFileIntoByteArray(file)).split("\n");
        for (final String line : evidenceTargetFileLines) {
            links.add(EvidenceTargetLink.fromBedpeString(line.trim(), dictionary));
        }
        return links;
    }

}
