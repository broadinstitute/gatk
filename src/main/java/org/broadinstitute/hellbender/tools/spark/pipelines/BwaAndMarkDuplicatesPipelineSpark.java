package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkPipelineProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

/**
 * Runs BWA and MarkDuplicates on Spark. It's an example of how to compose those two tools.
 */
@CommandLineProgramProperties(
        summary = "Takes name-sorted file and runs BWA and MarkDuplicates.",
        oneLineSummary = "Takes name-sorted file and runs BWA and MarkDuplicates.",
        usageExample = "BwaAndMarkDuplicatesPipelineSpark -I single.bam -R referenceURL -O file:///tmp/output.bam",
        programGroup = SparkPipelineProgramGroup.class
)
public final class BwaAndMarkDuplicatesPipelineSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    /**
     *  command-line arguments to specify the BWA behavior
     */
    @ArgumentCollection
    private BwaArgumentCollection bwaArgs = new BwaArgumentCollection();

    @Argument(shortName = "DS", fullName ="duplicates_scoring_strategy", doc = "The scoring strategy for choosing the non-duplicate among candidates.")
    public MarkDuplicatesScoringStrategy duplicatesScoringStrategy = MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String output;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> initialReads = getReads();
        final String referenceFileName = referenceArguments.getReferenceFileName();
        final BwaSparkEngine engine = new BwaSparkEngine(bwaArgs.numThreads, bwaArgs.fixedChunkSize, referenceFileName);
        final SAMFileHeader readsHeader = engine.makeHeaderForOutput(getHeaderForReads(), getReferenceSequenceDictionary());
        final JavaRDD<GATKRead> alignedReads = engine.alignWithBWA(ctx, initialReads, readsHeader);

        final JavaRDD<GATKRead> markedReadsWithOD = MarkDuplicatesSpark.mark(alignedReads, getHeaderForReads(), duplicatesScoringStrategy, new OpticalDuplicateFinder(), getRecommendedNumReducers());
        final JavaRDD<GATKRead> markedReads = MarkDuplicatesSpark.cleanupTemporaryAttributes(markedReadsWithOD);
        writeReads(ctx, output, markedReads);
    }
}
