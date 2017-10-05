package org.broadinstitute.hellbender.tools.spark.bwa;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;

import java.io.IOException;

@CommandLineProgramProperties(summary = "Runs BWA",
        oneLineSummary = "BWA on Spark",
        programGroup = SparkProgramGroup.class)
@BetaFeature
public final class BwaSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "the output bam",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String output;

    @ArgumentCollection
    public final BwaArgumentCollection bwaArgs = new BwaArgumentCollection();

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        try ( final BwaSparkEngine engine =
                      new BwaSparkEngine(ctx, referenceArguments.getReferenceFileName(), bwaArgs.indexImageFile, getHeaderForReads(), getReferenceSequenceDictionary()) ) {
            final JavaRDD<GATKRead> reads = !bwaArgs.singleEndAlignment ? engine.alignPaired(getReads()) : engine.alignUnpaired(getReads());

            try {
                ReadsSparkSink.writeReads(ctx, output, null, reads, engine.getHeader(),
                                            shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE);
            } catch (final IOException e) {
                throw new GATKException("Unable to write aligned reads", e);
            }
        }
    }
}
