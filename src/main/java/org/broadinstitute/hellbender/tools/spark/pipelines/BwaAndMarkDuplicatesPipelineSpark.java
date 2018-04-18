package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

import java.io.IOException;

/**
 * Runs BWA and MarkDuplicates on Spark. It's an example of how to compose those two tools.
 */
@CommandLineProgramProperties(
        summary = "Takes name-sorted file and runs BWA and MarkDuplicates.",
        oneLineSummary = "Takes name-sorted file and runs BWA and MarkDuplicates.",
        usageExample = "BwaAndMarkDuplicatesPipelineSpark -I single.bam -R referenceURL -O file:///tmp/output.bam",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class BwaAndMarkDuplicatesPipelineSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @ArgumentCollection
    public final BwaArgumentCollection bwaArgs = new BwaArgumentCollection();

    @Argument(shortName = "DS", fullName ="duplicates_scoring_strategy", doc = "The scoring strategy for choosing the non-duplicate among candidates.")
    public MarkDuplicatesScoringStrategy duplicatesScoringStrategy = MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    protected String output;

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        try (final BwaSparkEngine bwaEngine = new BwaSparkEngine(ctx, referenceArguments.getReferenceFileName(), bwaArgs.indexImageFile, getHeaderForReads(), getReferenceSequenceDictionary())) {
            final ReadFilter filter = makeReadFilter(bwaEngine.getHeader());
            final JavaRDD<GATKRead> alignedReads = bwaEngine.alignPaired(getUnfilteredReads()).filter(filter::test);
            final JavaRDD<GATKRead> markedReads = MarkDuplicatesSpark.mark(alignedReads, bwaEngine.getHeader(), duplicatesScoringStrategy, new OpticalDuplicateFinder(), getRecommendedNumReducers());
            try {
                ReadsSparkSink.writeReads(ctx, output,
                        referenceArguments.getReferencePath().toAbsolutePath().toUri().toString(),
                        markedReads, bwaEngine.getHeader(),
                        shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE,
                        getRecommendedNumReducers());
            } catch (IOException e) {
                throw new GATKException("unable to write bam: " + e);
            }
        }
    }
}
