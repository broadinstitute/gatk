package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.AddContextDataToReadSpark;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.JoinStrategy;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.transforms.BaseRecalibratorSparkFn;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;

import java.io.PrintStream;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).",
        oneLineSummary = "BaseRecalibrator on Spark",
        programGroup = SparkProgramGroup.class
)
@DocumentedFeature
public class BaseRecalibratorSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @Override
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return BaseRecalibrator.getStandardBQSRReadFilterList();
    }

    @Argument(doc = "the known variants", shortName = "knownSites", fullName = "knownSites", optional = false)
    private List<String> knownVariants;

    @Argument(doc = "the join strategy for reference bases and known variants", shortName = "joinStrategy", fullName = "joinStrategy", optional = true)
    private JoinStrategy joinStrategy = JoinStrategy.BROADCAST;

    @Argument(doc = "Path to save the final recalibration tables to.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputTablesPath = null;

    /**
     * all the command line arguments for BQSR and its covariates
     */
    @ArgumentCollection(doc = "all the command line arguments for BQSR and its covariates")
    private final RecalibrationArgumentCollection bqsrArgs = new RecalibrationArgumentCollection();

    @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases. Only applies when using the OVERLAPS_PARTITIONER join strategy.", optional = true)
    public int readShardSize = 10000;

    @Argument(fullName="readShardPadding", shortName="readShardPadding", doc = "Each read shard has this many bases of extra context on each side. Only applies when using the OVERLAPS_PARTITIONER join strategy.", optional = true)
    public int readShardPadding = 1000;

    @Override
    protected void runTool( JavaSparkContext ctx ) {
        if (joinStrategy == JoinStrategy.BROADCAST && ! getReference().isCompatibleWithSparkBroadcast()){
            throw new UserException.Require2BitReferenceForBroadcast();
        }

        JavaRDD<GATKRead> initialReads = getReads();

        if (joinStrategy.equals(JoinStrategy.OVERLAPS_PARTITIONER)) {
            // the overlaps partitioner requires that reads are coordinate-sorted
            final SAMFileHeader readsHeader = getHeaderForReads().clone();
            readsHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            initialReads = SparkUtils.coordinateSortReads(initialReads, readsHeader, numReducers);
        }

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        JavaRDD<GATKVariant> bqsrKnownVariants = variantsSparkSource.getParallelVariants(knownVariants, getIntervals());

        // TODO: Look into broadcasting the reference to all of the workers. This would make AddContextDataToReadSpark
        // TODO: and ApplyBQSRStub simpler (#855).
        JavaPairRDD<GATKRead, ReadContextData> rddReadContext = AddContextDataToReadSpark.add(ctx, initialReads, getReference(), bqsrKnownVariants, joinStrategy, getReferenceSequenceDictionary(), readShardSize, readShardPadding);

        // TODO: broadcast the reads header?
        final RecalibrationReport bqsrReport = BaseRecalibratorSparkFn.apply(rddReadContext, getHeaderForReads(), getReferenceSequenceDictionary(), bqsrArgs);

        try ( final PrintStream reportStream = new PrintStream(BucketUtils.createFile(outputTablesPath, getAuthenticatedGCSOptions())) ) {
            RecalUtils.outputRecalibrationReport(reportStream, bqsrArgs, bqsrReport.getQuantizationInfo(), bqsrReport.getRecalibrationTables(), bqsrReport.getCovariates());
        }
    }
}
