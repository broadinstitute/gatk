package org.broadinstitute.hellbender.tools.spark.bwa;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.SequenceDictionaryValidationArgumentCollection;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

@DocumentedFeature
@CommandLineProgramProperties(summary = "Align reads using BWA",
        oneLineSummary = "Align reads to a given reference using BWA on Spark",
        programGroup = ReadDataManipulationProgramGroup.class)
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
    public List<ReadFilter> getDefaultReadFilters() {
        // 1) unmapped or neither secondary nor supplementary and 2) has some sequence
        return Arrays.asList(ReadFilterLibrary.PRIMARY_LINE, ReadFilterLibrary.SEQ_IS_STORED);
    }

    @Override
    public SequenceDictionaryValidationArgumentCollection getSequenceDictionaryValidationArgumentCollection(){
        return new SequenceDictionaryValidationArgumentCollection.NoValidationCollection();
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        try ( final BwaSparkEngine bwaEngine =
                      new BwaSparkEngine(ctx, referenceArguments.getReferenceFileName(), bwaArgs.indexImageFile, getHeaderForReads(), getReferenceSequenceDictionary()) ) {
            final JavaRDD<GATKRead> reads;
            if (bwaArgs.singleEndAlignment) {
                reads = bwaEngine.alignUnpaired(getReads());
            } else {
                // filter reads after alignment in the case of paired reads since filtering does not know about pairs
                final ReadFilter filter = makeReadFilter(bwaEngine.getHeader());
                reads = bwaEngine.alignPaired(getUnfilteredReads()).filter(filter::test);
            }
            try {
                ReadsSparkSink.writeReads(ctx, output, null, reads, bwaEngine.getHeader(),
                                            shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE);
            } catch (final IOException e) {
                throw new GATKException("Unable to write aligned reads", e);
            }
        }
    }
}
