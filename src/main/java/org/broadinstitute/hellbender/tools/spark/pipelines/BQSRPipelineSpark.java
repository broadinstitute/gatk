package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.spark.JoinReadsWithVariants;
import org.broadinstitute.hellbender.tools.ApplyBQSRUniqueArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.transforms.ApplyBQSRSparkFn;
import org.broadinstitute.hellbender.tools.spark.transforms.BaseRecalibratorSparkFn;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.util.List;

/**
 * The full BQSR pipeline in one tool to run on Spark.
 * The final result is analysis-ready reads.
 * This runs BaseRecalibrator and then ApplyBQSR to give a BAM with recalibrated base qualities.
 *
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A BAM or CRAM file containing input read data</li>
 *     <li> A database of known polymorphic sites to skip over.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <p> A BAM or CRAM file containing the recalibrated read data</p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk BQSRPipelineSpark \
 *   -R gs://my-gcs-bucket/reference.fasta \
 *   -I gs://my-gcs-bucket/input.bam \
 *   --known-sites gs://my-gcs-bucket/sites_of_variation.vcf \
 *   --known-sites gs://my-gcs-bucket/another/optional/setOfSitesToMask.vcf \
 *   -O gs://my-gcs-bucket/output.bam \
 *   -- \
 *   --sparkRunner GCS \
 *   --cluster my-dataproc-cluster
 * </pre>
 */
@CommandLineProgramProperties(
        summary = BQSRPipelineSpark.USAGE_SUMMARY,
        oneLineSummary = BQSRPipelineSpark.USAGE_ONE_LINE_SUMMARY,
        usageExample = "BQSRPipelineSpark -I in.bam --known-sites in.vcf -O out.bam",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class BQSRPipelineSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    static final String USAGE_ONE_LINE_SUMMARY = "Both steps of BQSR (BaseRecalibrator and ApplyBQSR) on Spark";
    static final String USAGE_SUMMARY = "This tools performs 2 steps of BQSR - " +
            "creation of recalibration tables and rewriting of the bam, " +
            "without writing the tables to disk. ";

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @Argument(doc = "the known variants", fullName = BaseRecalibrator.KNOWN_SITES_ARG_FULL_NAME, optional = false)
    protected List<String> knownVariants;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String output;

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection bqsrArgs = new RecalibrationArgumentCollection();

    /**
     * command-line arguments to fine tune the apply BQSR step.
     */
    @ArgumentCollection
    public ApplyBQSRUniqueArgumentCollection applyBqsrArgs = new ApplyBQSRUniqueArgumentCollection();

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        String referenceFileName = addReferenceFilesForSpark(ctx, referenceArguments.getReferencePath());
        List<String> localKnownSitesFilePaths = addVCFsForSpark(ctx, knownVariants);

        //Should this get the getUnfilteredReads? getReads will merge default and command line filters.
        //but the code below uses other filters for other parts of the pipeline that do not honor
        //the commandline.
        final JavaRDD<GATKRead> initialReads = getReads();

        // The initial reads have already had the WellformedReadFilter applied to them, which
        // is all the filtering that ApplyBQSR wants. BQSR itself wants additional filtering
        // performed, so we do that here.
        //NOTE: this filter doesn't honor enabled/disabled commandline filters
        final ReadFilter bqsrReadFilter = ReadFilter.fromList(BaseRecalibrator.getBQSRSpecificReadFilterList(), getHeaderForReads());
        final JavaRDD<GATKRead> filteredReadsForBQSR = initialReads.filter(read -> bqsrReadFilter.test(read));

        JavaPairRDD<GATKRead, Iterable<GATKVariant>> readsWithVariants = JoinReadsWithVariants.join(filteredReadsForBQSR, localKnownSitesFilePaths);
        //note: we use the reference dictionary from the reads themselves.
        final RecalibrationReport bqsrReport = BaseRecalibratorSparkFn.apply(readsWithVariants, getHeaderForReads(), referenceFileName, bqsrArgs);

        final Broadcast<RecalibrationReport> reportBroadcast = ctx.broadcast(bqsrReport);
        final JavaRDD<GATKRead> finalReads = ApplyBQSRSparkFn.apply(initialReads, reportBroadcast, getHeaderForReads(), applyBqsrArgs.toApplyBQSRArgumentCollection(bqsrArgs));

        writeReads(ctx, output, finalReads);
    }
}
