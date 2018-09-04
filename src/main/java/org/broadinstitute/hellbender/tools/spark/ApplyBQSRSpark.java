package org.broadinstitute.hellbender.tools.spark;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.transforms.ApplyBQSRSparkFn;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

/**
 * Apply base quality score recalibration with Spark.
 * <p>
 * <p>This tool performs the second pass in a two-stage process called Base Quality Score Recalibration (BQSR).
 * Specifically, it recalibrates the base qualities of the input reads based on the recalibration table produced by
 * the BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file.</p>
 * <p>See <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.
 * </p>
 * <h3>Usage examples</h3>
 * <pre>
 * gatk ApplyBQSRSpark \
 * -I gs://my-gcs-bucket/input.bam \
 * -bqsr gs://my-gcs-bucket/recalibration.table \
 * -O gs://my-gcs-bucket/output.bam \
 * -- \
 * --sparkRunner GCS \
 * --cluster my-dataproc-cluster
 * </pre>
 * <p>
 * To additionally bin base qualities:
 * <pre>
 * gatk ApplyBQSRSpark \
 * -I gs://my-gcs-bucket/input.bam \
 * -bqsr gs://my-gcs-bucket/recalibration.table \
 * --static-quantized-quals 10 --static-quantized-quals 20 \
 * --static-quantized-quals 30 --static-quantized-quals 40 \
 * -O gs://my-gcs-bucket/output.bam \
 * -- \
 * --sparkRunner GCS \
 * --cluster my-dataproc-cluster
 * </pre>
 */

@CommandLineProgramProperties(
        summary=ApplyBQSRSpark.USAGE_SUMMARY,
        oneLineSummary=ApplyBQSRSpark.USAGE_ONE_LINE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class
)

@DocumentedFeature
@BetaFeature
public final class ApplyBQSRSpark extends GATKSparkTool {
    private static final long serialVersionUID = 0l;
    static final String USAGE_ONE_LINE_SUMMARY = "Apply base quality score recalibration on Spark";
    static final String USAGE_SUMMARY = "Apply a linear base quality recalibration model trained with the BaseRecalibrator tool on Spark.";

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
        Broadcast<RecalibrationReport> recalibrationReportBroadCast = ctx.broadcast(new RecalibrationReport(BucketUtils.openFile(bqsrRecalFile)));
        final JavaRDD<GATKRead> recalibratedReads = ApplyBQSRSparkFn.apply(initialReads, recalibrationReportBroadCast, getHeaderForReads(), applyBQSRArgs);
        writeReads(ctx, output, recalibratedReads);
    }
}
