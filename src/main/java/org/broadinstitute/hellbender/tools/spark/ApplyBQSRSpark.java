package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.tools.spark.transforms.ApplyBQSRSparkFn;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.util.List;

@CommandLineProgramProperties(summary="Apply Base Quality Recalibration to a bam file using spark",
        oneLineSummary="apply BQSR on spark",
        programGroup = SparkProgramGroup.class)
public final class ApplyBQSRSpark extends SparkCommandLineProgram{
    private static final long serialVersionUID = 0l;

    @Argument(doc = "uri for the input bam, either a local file path or a gs:// bucket path",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false)
    private String bam;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    /**
     * Enables recalibration of base qualities.
     * The covariates tables are produced by the BaseRecalibrator tool.
     * Please be aware that you should only run recalibration with the covariates file created on the same input bam(s).
     */
    @Argument(fullName="bqsr_recal_file", shortName="bqsr", doc="Input covariates table file for base quality score recalibration")
    private String bqsrRecalFile;

    @ArgumentCollection
    private IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @ArgumentCollection
    private ApplyBQSRArgumentCollection applyBQSRArgs = new ApplyBQSRArgumentCollection();

    @Override
    protected void runPipeline(JavaSparkContext ctx) {
        ReadsSparkSource readSource = new ReadsSparkSource(ctx);
        SAMFileHeader readsHeader = ReadsSparkSource.getHeader(ctx, bam);
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary())
                : IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());
        JavaRDD<GATKRead> initialReads = readSource.getParallelReads(bam, intervals);

        Broadcast<RecalibrationReport> recalibrationReportBroadCast = ctx.broadcast(new RecalibrationReport(BucketUtils.openFile(bqsrRecalFile, apiKey)));
        final JavaRDD<GATKRead> recalibratedReads = ApplyBQSRSparkFn.apply(initialReads, recalibrationReportBroadCast, readsHeader, applyBQSRArgs);

        try {
            ReadsSparkSink.writeReads(ctx, output, recalibratedReads, readsHeader, true);
        } catch (IOException e) {
            throw new GATKException("unable to write bam: " + e);
        }
    }

    @Override
    protected String getProgramName() {
        return "ApplyBQSRSpark";
    }
}
