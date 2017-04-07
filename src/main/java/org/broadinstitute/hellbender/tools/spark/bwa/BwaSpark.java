package org.broadinstitute.hellbender.tools.spark.bwa;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

@CommandLineProgramProperties(summary = "Runs BWA",
        oneLineSummary = "BWA on Spark",
        programGroup = SparkProgramGroup.class)
public final class BwaSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String output;

    @Argument(doc = "the bwa mem index image file name that you've distributed to each executor",
            fullName = "bwamemIndexImage")
    private String indexImageFile;

    @Argument(doc = "do single-ended alignment", fullName = "singleEndedMode", optional = true)
    private boolean singleEndedMode = false;

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
    protected void runTool( final JavaSparkContext ctx ) {
        try ( final BwaSparkEngine engine =
                      new BwaSparkEngine(ctx,
                              indexImageFile,
                              getHeaderForReads(),
                              getReferenceSequenceDictionary(),
                              !singleEndedMode) ) {
            final JavaRDD<GATKRead> reads = engine.align(getReads());

            try {
                ReadsSparkSink.writeReads(ctx, output, null, reads, engine.getHeader(),
                        shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE);
            } catch ( final IOException e ) {
                throw new GATKException("Unable to write aligned reads", e);
            }
        }
    }
}
