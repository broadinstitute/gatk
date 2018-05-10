package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.utils.FlatMapGluer;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import static org.broadinstitute.hellbender.tools.spark.sv.evidence.FindBreakpointEvidenceSpark.buildMetadata;
import static org.broadinstitute.hellbender.tools.spark.sv.evidence.FindBreakpointEvidenceSpark.findGenomewideHighCoverageIntervalsToIgnore;

/**
 * (Internal) Extracts evidence of structural variations from reads
 *
 * <p>This tool is used in development and should not be of interest to most researchers.  It repackages the first
 * two steps of the structural variation workflow as a separate tool for the convenience of developers.</p>
 * <p>This tool examines a SAM/BAM/CRAM for reads, or groups of reads, that demonstrate evidence of a structural
 * variation in the vicinity.  It records this evidence as a group of text files in a specified output directory
 * on Spark's HDFS file system.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>A file of paired-end, aligned and coordinate-sorted reads.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A directory of text files describing the evidence for structural variation discovered.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk ExtractSVEvidenceSpark \
 *     -I input_reads.bam \
 *     -O hdfs://my_cluster-m:8020/output_directory
 *     --aligner-index-image ignored --kmers-to-ignore ignored
 * </pre>
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Extracts evidence of structural variations from reads",
        summary =
        "This tool is used in development and should not be of interest to most researchers.  It packages one step" +
        " of the structural variation workflow as a separate tool for the convenience of developers." +
        " This tool examines a SAM/BAM/CRAM for reads, or groups of reads, that demonstrate evidence of a structural" +
        " variation in the vicinity.  It records this evidence as a group of text files in a specified output directory" +
        " on Spark's HDFS file system.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public final class ExtractSVEvidenceSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params =
            new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();


    @Argument(doc = "HDFS path for output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputDir;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        final JavaRDD<GATKRead> unfilteredReads = getUnfilteredReads();
        final SVReadFilter filter = new SVReadFilter(params);
        final ReadMetadata readMetadata = buildMetadata(params, getHeaderForReads(), unfilteredReads, filter, logger);
        final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadata);
        final int allowedOverhang = params.allowedShortFragmentOverhang;
        final int minEvidenceMapQ = params.minEvidenceMapQ;

        final SVIntervalTree<SVInterval> highCoverageSubintervalTree = findGenomewideHighCoverageIntervalsToIgnore(params,
                readMetadata, ctx, getHeaderForReads(), unfilteredReads, filter, logger, broadcastMetadata);
        final Broadcast<SVIntervalTree<SVInterval>> broadcastHighCoverageSubIntervals = ctx.broadcast(highCoverageSubintervalTree);

        unfilteredReads
            .mapPartitions(readItr -> {
                final GATKRead sentinel = new SAMRecordToGATKReadAdapter(null);
                return FlatMapGluer.applyMapFunc(
                        new ReadClassifier(broadcastMetadata.value(),
                                sentinel,
                                allowedOverhang,
                                filter,
                                broadcastHighCoverageSubIntervals.getValue()),
                        readItr, sentinel);
            }, true)
            .map(e -> e.stringRep(broadcastMetadata.getValue(), minEvidenceMapQ))
            .saveAsTextFile(outputDir);
    }
}
