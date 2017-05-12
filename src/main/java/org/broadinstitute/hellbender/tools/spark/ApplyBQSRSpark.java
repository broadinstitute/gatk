package org.broadinstitute.hellbender.tools.spark;

import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.transforms.ApplyBQSRSparkFn;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;

@CommandLineProgramProperties(summary="Apply Base Quality Recalibration to a SAM/BAM/CRAM file using Spark",
        oneLineSummary="ApplyBQSR on Spark",
        programGroup = SparkProgramGroup.class)
@DocumentedFeature
public final class ApplyBQSRSpark extends GATKSparkTool {
    private static final long serialVersionUID = 0l;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    /**
     * Enables recalibration of base qualities.
     * The covariates tables are produced by the BaseRecalibrator tool.
     * Please be aware that you should only run recalibration with the covariates file created on the same input bam(s).
     */
    @Argument(fullName= StandardArgumentDefinitions.BQSR_TABLE_LONG_NAME, shortName= StandardArgumentDefinitions.BQSR_TABLE_SHORT_NAME, doc="Input covariates table file for base quality score recalibration")
    private String bqsrRecalFile;

    @ArgumentCollection
    private ApplyBQSRArgumentCollection applyBQSRArgs = new ApplyBQSRArgumentCollection();

    @Override
    protected void runTool(JavaSparkContext ctx) {
        JavaRDD<GATKRead> initialReads = getReads();
        final GCSOptions gcsOptions = getAuthenticatedGCSOptions(); // null if we have no api key
        Broadcast<RecalibrationReport> recalibrationReportBroadCast = ctx.broadcast(new RecalibrationReport(BucketUtils.openFile(bqsrRecalFile)));
        final JavaRDD<GATKRead> recalibratedReads = ApplyBQSRSparkFn.apply(initialReads, recalibrationReportBroadCast, getHeaderForReads(), applyBQSRArgs);
        writeReads(ctx, output, recalibratedReads);
    }
}
