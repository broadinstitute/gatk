package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.PathSeqProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import scala.Tuple2;

import java.io.IOException;

/**
 * This Spark tool is the first step in the PathSeq pipeline.
 * Input: set of unaligned or host-aligned reads. Optional: host K-mer file and bwa mem image.
 * Output: set of high quality non-host reads.
 *
 * Filtering steps:
 * 1)  Remove secondary and supplementary alignments
 * 2)  Trim adapter sequences
 * 3)  Mask sequences with excessive A/T or G/C content
 * 4)  Mask repetitive sequences with 'N' and base quality --dustPhred using symmetric DUST
 * 5)  Hard clip read ends using base qualities
 * 6)  Remove reads shorter than --minClippedReadLength
 * 7)  Mask bases whose Phred score is less than --minBaseQuality with 'N'
 * 8)  Remove reads whose fraction of bases that are 'N' is greater than --maxAmbiguousBases
 * 9)  If specified, remove reads containing one or more kmers --kmerLibraryPath
 * 10) If specified, remove reads that align to the host BWA image --filterBwaImage with at least --minCoverage and --minIdentity
 * 11) If --filterDuplicates is set, remove exact duplicate reads (not using Mark Duplicates because it requires aligned reads)
 *
 * Notes:
 *
 * - Steps 2 - 6 can be skipped by setting --skipFilters.
 * - The tool assumes the BAM file is unaligned by default. If the BAM is aligned to the host, use --isHostAligned to filter.
 * - Output will be two BAM files, outputPath.paired.bam and outputPath.unpaired.bam containing paired and unpaired reads.
 * - If the resulting set of reads is empty, the file will not be written.
 * - If --unpaired is set, pairedness flags will not be corrected after filtering, and all reads will be written to
 *     outputPath.unpaired.bam. This improves performance by avoiding shuffles but violates the SAM format specification.
 *
 */
@CommandLineProgramProperties(summary = "First, low-quality and repetitive sequences reads are filtered: \n\t(1) Remove secondary and " +
        "supplementary alignments; \n\t(2) trim adapter sequences; \n\t(3) mask sequences with excessive A/T or G/C content; " +
        "\n\t(4) mask repetitive sequences using the sDUST algorithm; \n\t(5) hard clip according to masked " +
        "and low-quality bases; \n\t(6) remove reads below minimum length; \n\t(7) mask low-quality bases; \n\t(8) filter reads with" +
        " too many masked bases. \nHost reads are then filtered using optionally-supplied host k-mer database (created " +
        "with PathSeqBuildKmers) and host BWA index image (created with BwaMemIndexImageCreator). Lastly, exact " +
        "duplicate sequences are removed.",
        oneLineSummary = "Step 1: Filters low-quality, low-complexity, duplicate, and host reads",
        programGroup = PathSeqProgramGroup.class)
@BetaFeature
public final class PathSeqFilterSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "Base uri for the output file(s). Paired and unpaired reads will be written to uri appended with" +
            " '.paired.bam' and '.unpaired.bam'",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public String outputPath;

    @ArgumentCollection
    public PSFilterArgumentCollection filterArgs = new PSFilterArgumentCollection();

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        filterArgs.doReadFilterArgumentWarnings(getCommandLineParser().getPluginDescriptor(GATKReadFilterPluginDescriptor.class), logger);
        final SAMFileHeader header = PSUtils.checkAndClearHeaderSequences(getHeaderForReads(), filterArgs, logger);

        final PSFilter filter = new PSFilter(ctx, filterArgs, getReads(), header);

        final Tuple2<JavaRDD<GATKRead>, JavaRDD<GATKRead>> result = filter.doFilter();
        final JavaRDD<GATKRead> pairedReads = result._1;
        final JavaRDD<GATKRead> unpairedReads = result._2;

        if (!pairedReads.isEmpty()) {
            header.setSortOrder(SAMFileHeader.SortOrder.queryname);
            writeReads(ctx, outputPath + ".paired.bam", pairedReads, header);
        } else {
            logger.info("No paired reads to write - BAM will not be written.");
        }
        if (!unpairedReads.isEmpty()) {
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
            writeReads(ctx, outputPath + ".unpaired.bam", unpairedReads, header);
        } else {
            logger.info("No unpaired reads to write - BAM will not be written.");
        }
        filter.close();
    }

}
