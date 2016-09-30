package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.Collections;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary="Compute metrics on the number of chimeric read pairs that cluster together",
        oneLineSummary="Compute metrics on the number of chimeric read pairs that cluster together",
        programGroup = SparkProgramGroup.class)
public class CalculateChimericFreqRateSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final SAMSequenceDictionary referenceSequenceDictionary = getReferenceSequenceDictionary();
        final JavaPairRDD<Tuple2<Integer, Integer>, Integer> binPairCounts = getReads().filter(r -> !r.isSecondaryAlignment() && !r.isSupplementaryAlignment() && !r.failsVendorQualityCheck() && r.isFirstOfPair() && !r.isUnmapped() && !r.mateIsUnmapped()).
                map(r -> getLinkedBins(r, referenceSequenceDictionary)).
                filter(p -> Math.abs(p._1() - p._2()) > 1).mapToPair(p -> new Tuple2<>(p, 1)).reduceByKey((x, y) -> x + y);
        final Map<Integer, Long> freqMap = binPairCounts.map(Tuple2::_2).countByValue();

        freqMap.entrySet().forEach(v -> System.out.println(v.getKey() + "\t" + v.getValue()));
    }



    private static Tuple2<Integer, Integer> getLinkedBins(final GATKRead r, final SAMSequenceDictionary referenceSequenceDictionary) {
        final String contig = r.getContig();
        final int contigIndex = referenceSequenceDictionary.getSequenceIndex(contig);
        final int mateContigIndex = referenceSequenceDictionary.getSequenceIndex(r.getMateContig());

        return new Tuple2<>(contigIndex * 10000000 + r.getStart() / 500,  mateContigIndex * 10000000 + r.getMateStart() / 500);
    }
}
