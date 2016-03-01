package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import scala.Tuple2;

import java.io.IOException;
import java.util.ArrayList;

import java.util.List;

@CommandLineProgramProperties(summary="Gather clustered split reads using spark, simple version",
        oneLineSummary="Gather clustered split reads using sparks, simple version",
        programGroup = SparkProgramGroup.class)
public class GatherSplitReadsSimpleSpark extends GATKSparkTool
{
    private static final long serialVersionUID = 1L;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public ReadFilter makeReadFilter() {
        return ReadFilterLibrary.ALLOW_ALL_READS;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> reads = getReads();

        System.err.println("about to run");

        /**
         * Generate mapping between (chromosome + left & right soft-clipping positions, the read) in the read (two entries of both exist),
         *     then group by (chromosome + left & right soft-clipping positions)
         *     then filter out if less than 2 reads support this particular clipping location
         *     then return all supportive reads of filter-passing clipping locations.
         */
        final JavaRDD<GATKRead> clusteredSplitReads = reads.flatMapToPair(new SRPairFlatMapFunction())
                                                           .groupByKey()
                                                           .filter(tuple -> { final Iterable<GATKRead> values = tuple._2();
                                                                              int valueCount = 0;
                                                                              for (GATKRead value : values) {
                                                                                  valueCount++;
                                                                                  if (valueCount >= 2) {
                                                                                      return true;
                                                                                  }
                                                                              }
                                                                              return false;
                                                                              })
                                                           .flatMap(Tuple2::_2);

        try {
            ReadsSparkSink.writeReads(ctx, output, clusteredSplitReads, getHeaderForReads(), ReadsWriteFormat.SINGLE);
        } catch (IOException e) {
            throw new GATKException("unable to write bam: " + e);
        }

    }

    /**
     * Cigar-string based method for mapping (chromosome + left & right soft-clipping positions) to read.
     * @param read
     * @param out
     */
    private static void addReadKeyedByClippingLocationsToOut(final GATKRead read, final List<Tuple2<SimpleInterval, GATKRead>> out) {
        final Cigar cigar = read.getCigar();
        // first CigarElement is the first non-hard-clipped element
        final int firstCigarElement = cigar.getCigarElement(0).getOperator() == CigarOperator.HARD_CLIP ? 1 : 0;
        final int firstElementLength = cigar.getCigarElement(firstCigarElement).getLength();
        if (cigar.getCigarElement(firstCigarElement).getOperator() == CigarOperator.SOFT_CLIP &&
                firstElementLength >= 40) {
            final SimpleInterval leftClipLoc = new SimpleInterval(read.getContig(), read.getStart(), read.getStart());
            out.add(new Tuple2<>(leftClipLoc, read));
        }
        // similarly for right-end
        final int lastCigarElement = cigar.getCigarElement(cigar.numCigarElements() - 1).getOperator() == CigarOperator.HARD_CLIP ? cigar.numCigarElements() - 2 : cigar.numCigarElements() - 1;
        final int lastElementLength = cigar.getCigarElement(lastCigarElement).getLength();

        if (cigar.getCigarElement(lastCigarElement).getOperator() == CigarOperator.SOFT_CLIP &&
                lastElementLength >= 40) {
            final SimpleInterval leftClipLoc = new SimpleInterval(read.getContig(), read.getEnd(), read.getEnd());
            out.add(new Tuple2<>(leftClipLoc, read));
        }
    }

    private static boolean shouldConsiderRead(final GATKRead gatkRead) {
        return ! (gatkRead.failsVendorQualityCheck() || gatkRead.isUnmapped() || gatkRead.getMappingQuality() == 0);
    }

    private static class SRPairFlatMapFunction implements PairFlatMapFunction<GATKRead, SimpleInterval, GATKRead> {
        public static final long serialVersionUID = 1L;

        @Override
        public Iterable<Tuple2<SimpleInterval, GATKRead>> call(final GATKRead read) throws Exception {
            final List<Tuple2<SimpleInterval, GATKRead>> out = new ArrayList<>();
            if (shouldConsiderRead(read)) {
                addReadKeyedByClippingLocationsToOut(read, out);
            }
            return out;
        }
    }
}
