package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.read.markduplicates.DuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

import java.io.IOException;

@CommandLineProgramProperties(
        summary ="Marks duplicates on spark",
        oneLineSummary ="Mark Duplicates",
        programGroup = SparkProgramGroup.class)
public final class MarkDuplicatesSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String output;

    @Argument(doc = "Path to write duplication metrics to.", optional=true,
            shortName = "M", fullName = "METRICS_FILE")
    protected String metricsFile;

    @ArgumentCollection
    protected OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();

    @Argument(doc="The output parallelism, sets the number of reducers. Defaults to the number of partitions in the input.",
            shortName = "P", fullName = "parallelism", optional = true)
    protected int parallelism = 0;

    public static JavaRDD<GATKRead> mark(final JavaRDD<GATKRead> reads, final SAMFileHeader header,
                                         final OpticalDuplicateFinder opticalDuplicateFinder) {
        return mark(reads, header, opticalDuplicateFinder, reads.partitions().size());
    }

    public static JavaRDD<GATKRead> mark(final JavaRDD<GATKRead> reads, final SAMFileHeader header,
                                         final OpticalDuplicateFinder opticalDuplicateFinder, final int parallelism) {

        JavaRDD<GATKRead> primaryReads = reads.filter(v1 -> !isNonPrimary(v1));
        JavaRDD<GATKRead> nonPrimaryReads = reads.filter(v1 -> isNonPrimary(v1));
        JavaRDD<GATKRead> primaryReadsTransformed = MarkDuplicatesSparkUtils.transformReads(header, opticalDuplicateFinder, primaryReads, parallelism);

        return primaryReadsTransformed.union(nonPrimaryReads);
    }

    private static boolean isNonPrimary(GATKRead read) {
        return read.isSecondaryAlignment() || read.isSupplementaryAlignment() || read.isUnmapped();
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        JavaRDD<GATKRead> reads = getReads();
        if (parallelism == 0) { // use the number of partitions in the input
            parallelism = reads.partitions().size();
        }
        final OpticalDuplicateFinder finder = opticalDuplicatesArgumentCollection.READ_NAME_REGEX != null ?
                new OpticalDuplicateFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX, opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE, null) : null;

        final JavaRDD<GATKRead> finalReads = mark(reads, getHeaderForReads(), finder, parallelism);

        try {
            ReadsSparkSink.writeReads(ctx, output, finalReads, getHeaderForReads(), parallelism == 1 ? ReadsWriteFormat.SINGLE : ReadsWriteFormat.SHARDED);
        } catch (IOException e) {
            throw new GATKException("unable to write bam: " + e);
        }

        if (metricsFile != null) {
            final JavaPairRDD<String, DuplicationMetrics> metrics = MarkDuplicatesSparkUtils.generateMetrics(getHeaderForReads(), finalReads);
            MarkDuplicatesSparkUtils.saveMetricsRDD(metrics, metricsFile, getAuthHolder());
        }
    }
}
