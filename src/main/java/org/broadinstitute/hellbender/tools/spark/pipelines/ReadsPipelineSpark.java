package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
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
import org.broadinstitute.hellbender.tools.spark.bwa.BwaArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.tools.spark.transforms.ApplyBQSRSparkFn;
import org.broadinstitute.hellbender.tools.spark.transforms.BaseRecalibratorSparkFn;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;

import java.util.List;


@CommandLineProgramProperties(
        summary = "Takes aligned reads (likely from BWA) and runs MarkDuplicates and BQSR. The final result is analysis-ready reads.",
        oneLineSummary = "Takes aligned reads (likely from BWA) and runs MarkDuplicates and BQSR. The final result is analysis-ready reads",
        usageExample = "ReadsPipelineSpark -I single.bam -R referenceURL -knownSites variants.vcf -O file:///tmp/output.bam",
        programGroup = SparkPipelineProgramGroup.class
)

/**
 * ReadsPipelineSpark is our standard pipeline that takes aligned reads (likely from BWA) and runs MarkDuplicates
 * and BQSR. The final result is analysis-ready reads.
 */
public class ReadsPipelineSpark extends GATKSparkTool {
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

    @Argument(shortName = "DS", fullName ="duplicates_scoring_strategy", doc = "The scoring strategy for choosing the non-duplicate among candidates.")
    public MarkDuplicatesScoringStrategy duplicatesScoringStrategy = MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES;

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

    @Argument(shortName = "bwa", fullName = "run_bwa", doc = "run bwa alignment at the beginning of the pipeline", optional = true)
    public boolean runBWA = false;

    //HACK: this is only needed because BWA cannot use 2bit and this pipeline wants a 2bit reference
    @Argument(shortName = "bwaRef", fullName = "bwa_reference", doc = "Reference used by bwa", optional = true)
    public String referenceFileName;

    /**
     *  command-line arguments to specify the BWA behavior
     */
    @ArgumentCollection
    private BwaArgumentCollection bwaArgs = new BwaArgumentCollection();

    @Override
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        if (joinStrategy == JoinStrategy.BROADCAST && ! getReference().isCompatibleWithSparkBroadcast()){
            throw new UserException.Require2BitReferenceForBroadcast();
        }

        final SAMFileHeader header;
        final JavaRDD<GATKRead> alignedReads;
        if (runBWA) {
            if (referenceFileName == null){
                throw new UserException.MissingArgument("bwa_reference", "Provide a fasta file for bwa");
            }
            final JavaRDD<GATKRead> initialReads = getReads();

            //HACK: referenceFileName should be referenceArguments.getReferenceFileName(); once the 2bit issue is resolved
            final BwaSparkEngine bwaEngine = new BwaSparkEngine(bwaArgs.numThreads, bwaArgs.fixedChunkSize, referenceFileName);
            final SAMFileHeader readsHeader = bwaEngine.makeHeaderForOutput(getHeaderForReads(), getReferenceSequenceDictionary());
            alignedReads = bwaEngine.alignWithBWA(ctx, initialReads, readsHeader);
            header = readsHeader;
        } else {
            alignedReads = getReads();
            header = getHeaderForReads();
        }

        final JavaRDD<GATKRead> markedReadsWithOD = MarkDuplicatesSpark.mark(alignedReads, header, duplicatesScoringStrategy, new OpticalDuplicateFinder(), getRecommendedNumReducers());
        final JavaRDD<GATKRead> markedReads = MarkDuplicatesSpark.cleanupTemporaryAttributes(markedReadsWithOD);

        // The markedReads have already had the WellformedReadFilter applied to them, which
        // is all the filtering that MarkDupes and ApplyBQSR want. BQSR itself wants additional
        // filtering performed, so we do that here.
        final ReadFilter bqsrReadFilter = BaseRecalibrator.makeBQSRSpecificReadFilters();
        final JavaRDD<GATKRead> markedFilteredReadsForBQSR = markedReads.filter(bqsrReadFilter::test);

        final VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        final JavaRDD<GATKVariant> bqsrKnownVariants = variantsSparkSource.getParallelVariants(baseRecalibrationKnownVariants, getIntervals());

        final JavaPairRDD<GATKRead, ReadContextData> rddReadContext = AddContextDataToReadSpark.add(markedFilteredReadsForBQSR, getReference(), bqsrKnownVariants, joinStrategy);
        final RecalibrationReport bqsrReport = BaseRecalibratorSparkFn.apply(rddReadContext, header, getReferenceSequenceDictionary(), bqsrArgs);

        final Broadcast<RecalibrationReport> reportBroadcast = ctx.broadcast(bqsrReport);
        final JavaRDD<GATKRead> finalReads = ApplyBQSRSparkFn.apply(markedReads, reportBroadcast, header, applyBqsrArgs.toApplyBQSRArgumentCollection(bqsrArgs.PRESERVE_QSCORES_LESS_THAN));

        writeReads(ctx, output, finalReads);
    }
}
