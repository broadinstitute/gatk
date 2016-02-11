package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import scala.Tuple2;

import java.io.IOException;

@CommandLineProgramProperties(summary="Gather one-end-anchored alignments using spark",
        oneLineSummary="Gather one-end-anchored alignments using spark",
        programGroup = SparkProgramGroup.class)
public class GatherOneEndAnchoredPairsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final SAMFileHeader header = getHeaderForReads();

        final JavaRDD<GATKRead> reads = getReads();
        final JavaPairRDD<String, Iterable<GATKRead>> oeaReadsGroupedByName =
                reads
                        .filter(read ->
                                !read.failsVendorQualityCheck() &&
                                        ! read.isDuplicate() &&
                                        read.getMappingQuality() > 0 &&
                                        ! read.isSecondaryAlignment() &&
                                        ! read.isSupplementaryAlignment() &&
                                        (!read.isUnmapped() && read.mateIsUnmapped()) || (read.isUnmapped() && !read.mateIsUnmapped()))
                        .mapToPair(read -> new Tuple2<>(read.getName(), read))
                        .groupByKey();

        final JavaRDD<GATKRead> oeaReads = oeaReadsGroupedByName.flatMap(Tuple2::_2);

        try {
            ReadsSparkSink.writeReads(ctx,output,oeaReads,header, ReadsWriteFormat.SINGLE); }
        catch ( IOException e ) {
            throw new GATKException("Unable to write BAM" + output, e);
        }
    }

}
