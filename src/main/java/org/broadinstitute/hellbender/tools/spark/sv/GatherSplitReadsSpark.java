package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import scala.Tuple2;


import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(summary="Gather clustered split reads using spark",
        oneLineSummary="Gather clustered split reads using spark",
        programGroup = SparkProgramGroup.class)
public class GatherSplitReadsSpark extends GATKSparkTool {

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Override
    public boolean requiresReads() {
        return super.requiresReads();
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();

        // todo: group by key is really slow on large data, figure out another way
        final JavaRDD<GATKRead> clusteredSplitReads =
                reads.flatMapToPair(read -> {
                    final List<Tuple2<GenomeLoc, GATKRead>> out = new ArrayList<>();
                    if (shouldConsiderRead(read)) {
                        final Cigar cigar = read.getCigar();
                        final GenomeLocParser genomeLocParser = new GenomeLocParser(getHeaderForReads().getSequenceDictionary());
                        final int firstCigarElement = cigar.getCigarElement(0).getOperator() == CigarOperator.HARD_CLIP ? 1 : 0;
                        if (cigar.getCigarElement(firstCigarElement).getOperator() == CigarOperator.SOFT_CLIP &&
                                cigar.getCigarElement(firstCigarElement).getLength() >= 30) {
                            final GenomeLoc leftClipLoc = genomeLocParser.createGenomeLoc(read.getContig(), read.getStart(), read.getStart());
                            out.add(new Tuple2<>(leftClipLoc, read));
                        }
                        final int lastCigarElement = cigar.getCigarElement(cigar.numCigarElements() - 1).getOperator() == CigarOperator.HARD_CLIP ? cigar.numCigarElements() - 2 : cigar.numCigarElements() - 1;
                        if (cigar.getCigarElement(lastCigarElement).getOperator() == CigarOperator.SOFT_CLIP &&
                                cigar.getCigarElement(firstCigarElement).getLength() >= 30) {
                            final GenomeLoc leftClipLoc = genomeLocParser.createGenomeLoc(read.getContig(), read.getEnd(), read.getEnd());
                            out.add(new Tuple2<>(leftClipLoc, read));
                        }

                    }
                    return out;
                }).groupByKey().filter(tuple -> {
            final Iterable<GATKRead> values = tuple._2();
            int valueCount = 0;
            for (GATKRead value : values) {
                valueCount++;
                if (valueCount >= 2) {
                    return true;
                }
            }
            return false;
        }).flatMap(Tuple2::_2);

        try {
            ReadsSparkSink.writeReads(ctx, output, clusteredSplitReads, getHeaderForReads(), ReadsWriteFormat.SINGLE);
        } catch (IOException e) {
            throw new GATKException("unable to write bam: " + e);
        }

    }

    private boolean shouldConsiderRead(final GATKRead gatkRead) {
        return ! (gatkRead.failsVendorQualityCheck() || gatkRead.isUnmapped() || gatkRead.getMappingQuality() == 0);
    }
}
