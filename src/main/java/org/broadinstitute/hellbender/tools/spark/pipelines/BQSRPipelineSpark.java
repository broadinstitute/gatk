package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkPipelineProgramGroup;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.AddContextDataToReadSpark;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.JoinStrategy;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSRUniqueArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.transforms.ApplyBQSRSparkFn;
import org.broadinstitute.hellbender.tools.spark.transforms.BaseRecalibratorSparkFn;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;

import java.util.List;
@CommandLineProgramProperties(
        summary = "This tools performs 2 steps of BQSR - creation of recalibration tables and rewriting of the bam, without writing the tables to disk. ",
        oneLineSummary = "Both steps of BQSR (BaseRecalibrator and ApplyBQSR) on Spark",
        usageExample = "BQSRPipelineSpark -I in.bam --knownSites in.vcf -O out.bam",
        programGroup = SparkPipelineProgramGroup.class
)
/**
 * BQSR. The final result is analysis-ready reads.
 */
public final class BQSRPipelineSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @Argument(doc = "the known variants", shortName = "knownSites", fullName = "knownSites", optional = false)
    protected List<String> baseRecalibrationKnownVariants;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String output;

    @Argument(doc = "the join strategy for reference bases and known variants", shortName = "joinStrategy", fullName = "joinStrategy", optional = true)
    private JoinStrategy joinStrategy = JoinStrategy.BROADCAST;

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
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        if (joinStrategy == JoinStrategy.BROADCAST && ! getReference().isCompatibleWithSparkBroadcast()){
            throw new UserException.Require2BitReferenceForBroadcast();
        }
        //TOOO: should this get the getUnfilteredReads? getReads will apply default and command line filters
        final JavaRDD<GATKRead> initialReads = getReads();

        // The initial reads have already had the WellformedReadFilter applied to them, which
        // is all the filtering that ApplyBQSR wants. BQSR itself wants additional filtering
        // performed, so we do that here.
        //NOTE: this doesn't honor enabled/disabled commandline filters
        final ReadFilter bqsrReadFilter = BaseRecalibrator.makeBQSRSpecificReadFilters().stream()
                .map(f -> f.getFilterInstance(getHeaderForReads(), null))
                .reduce(ReadFilterLibrary.ALLOW_ALL_READS, (f1, f2) -> f1.and(f2));
        final JavaRDD<GATKRead> filteredReadsForBQSR = initialReads.filter(read -> bqsrReadFilter.test(read));

        final VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        final JavaRDD<GATKVariant> bqsrKnownVariants = variantsSparkSource.getParallelVariants(baseRecalibrationKnownVariants);

        final JavaPairRDD<GATKRead, ReadContextData> rddReadContext = AddContextDataToReadSpark.add(filteredReadsForBQSR, getReference(), bqsrKnownVariants, joinStrategy);
        //note: we use the reference dictionary from the reads themselves.
        final RecalibrationReport bqsrReport = BaseRecalibratorSparkFn.apply(rddReadContext, getHeaderForReads(), getHeaderForReads().getSequenceDictionary(), bqsrArgs);

        final Broadcast<RecalibrationReport> reportBroadcast = ctx.broadcast(bqsrReport);
        final JavaRDD<GATKRead> finalReads = ApplyBQSRSparkFn.apply(initialReads, reportBroadcast, getHeaderForReads(), applyBqsrArgs.toApplyBQSRArgumentCollection(bqsrArgs.PRESERVE_QSCORES_LESS_THAN));

        writeReads(ctx, output, finalReads);
    }
}
