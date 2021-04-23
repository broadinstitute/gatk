package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.List;


/**
 * SortSam on Spark (works on SAM/BAM/CRAM)
 *
 * <p>A Spark implementation of <a href='https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_SortSam.php'>Picard SortSam</a>. The Spark version can run in parallel on multiple cores on a local machine or multiple machines on a Spark cluster while still matching the output of the single-core Picard version. See <a href="https://software.broadinstitute.org/gatk/blog?id=23420">Blog#23420</a> for performance benchmarks.</p>
 *
 * <p>The tool sorts reads by coordinate order by default or alternatively by read name, the QNAME field, if asked with the '-SO queryname' option. The contig ordering in the reference dictionary defines coordinate order, and the tool uses the sequence dictionary represented by the @SQ header lines or that of the optionally provided reference to sort reads by the RNAME field. For those reads mapping to a contig, coordinate sorting further orders reads by the POS field of the SAM record, which contains the leftmost mapping position.</p>
 *
 *  <p>To queryname-sort, the tool first groups by readname and then deterministically sorts within a readname set by orientation, secondary and supplementary SAM flags. For paired-end reads, reads in the pair share the same queryname. Because aligners can generate secondary and supplementary alignments, queryname groups can consists of, e.g. more than two records for a paired-end pair.</p>
 *
 * <h3>Usage examples</h3>
 * Coordinate-sort aligned reads using all cores available locally
 * <pre>
 * gatk SortSamSpark \
 * -I aligned.bam \
 * -O coordinatesorted.bam
 * </pre>
 *
 * Queryname-sort reads using four cores on a Spark cluster
 * <pre>
 * gatk SortSamSpark \
 * -I coordinatesorted.bam \
 * -SO queryname \
 * -O querygroupsorted.bam \
 * -- \
 *  --spark-runner SPARK \
 *  --spark-master <SPARK-CLUSTER-NAME>\
 *  --num-executors 5 \
 *  --executor-cores 4
 * </pre>
 *
 * <h3>Notes</h3>
 * <ol>
 *     <li>This Spark tool requires a significant amount of disk operations. Run with both the input data and outputs on high throughput SSDs when possible. When pipelining this tool on Google Compute Engine instances, for best performance requisition machines with LOCAL SSDs.  </li>
 *     <li>Furthermore, we recommend explicitly setting the Spark temp directory to an available SSD when running this in local mode by adding the argument --conf 'spark.local.dir=/PATH/TO/TEMP/DIR'. See the discussion at https://gatkforums.broadinstitute.org/gatk/discussion/comment/56337 for details.</li>
 * </ol>
 */
@DocumentedFeature
@CommandLineProgramProperties(summary = "Sorts the input SAM/BAM/CRAM",
        oneLineSummary = "SortSam on Spark (works on SAM/BAM/CRAM)",
        programGroup = ReadDataManipulationProgramGroup.class)
@BetaFeature
public final class SortSamSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the output file path", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputFile;

    @Argument(doc="sort order of the output file", shortName = StandardArgumentDefinitions.SORT_ORDER_SHORT_NAME, fullName = StandardArgumentDefinitions.SORT_ORDER_LONG_NAME, optional = true)
    private SparkSortOrder sortOrder = SparkSortOrder.coordinate;

    /**
     * SortOrders that have corresponding implementations for spark.
     * These correspond to a subset of {@link SAMFileHeader.SortOrder}.
     */
    private enum SparkSortOrder {
        coordinate(SAMFileHeader.SortOrder.coordinate),
        queryname(SAMFileHeader.SortOrder.queryname);

        private final SAMFileHeader.SortOrder order;

        SparkSortOrder(SAMFileHeader.SortOrder order) {
            this.order = order;
        }

        public SAMFileHeader.SortOrder getSamOrder() {
             return order;
        }
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    protected void onStartup() {
        super.onStartup();
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();
        final int numReducers = getRecommendedNumReducers();
        logger.info("Using {} reducers", numReducers);

        final SAMFileHeader header = getHeaderForReads();
        header.setSortOrder(sortOrder.getSamOrder());

        writeReads(ctx, outputFile, reads, header, true);
    }
}
