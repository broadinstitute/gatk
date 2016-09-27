package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.FlatMapFunction2;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.iterators.IntervalOverlappingIterator;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Print read alignments in Pileup-style format on Spark",
        oneLineSummary = "Print read alignments in Pileup-style format on Spark",
        programGroup = SparkProgramGroup.class)
public final class PileupSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the output file path", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    @Argument(doc = "size to break intervals up into", shortName = "intervalShardSize", fullName = "intervalShardSize", optional = true)
    private int intervalShardSize = 10000;

    @Override
    public ReadFilter makeReadFilter(){
        return new WellformedReadFilter(getHeaderForReads())
                .and(ReadFilterLibrary.MAPPED)
                .and(ReadFilterLibrary.NOT_DUPLICATE)
                .and(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK)
                .and(ReadFilterLibrary.PRIMARY_ALIGNMENT);
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();

        final List<SimpleInterval> intervals = hasIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(getHeaderForReads().getSequenceDictionary());

        final List<SimpleInterval> intervalShards = getIntervalShards(intervals, intervalShardSize);

        JavaRDD<AlignmentContext> alignmentContexts;
        if (shuffle) {
            alignmentContexts = SparkUtils.joinOverlappingShuffle(ctx, reads, GATKRead.class, intervalShards,
                    new ShufflePileupFunction(getHeaderForReads()));
        } else {
            alignmentContexts = SparkUtils.joinOverlapping(ctx, reads, GATKRead.class, getBestAvailableSequenceDictionary(),
                    intervalShards, new PileupFunction(getHeaderForReads()));
        }

        Broadcast<ReferenceMultiSource> bReferenceSource = hasReference() ? ctx.broadcast(getReference()) : null;
        alignmentContexts.map(alignmentContext -> {
            // TODO: this reference lookup is potentially slow, push down to do once per interval shard
            byte[] bases = hasReference() ? bReferenceSource.getValue().getReferenceBases(null, new SimpleInterval(alignmentContext.getLocation())).getBases() : null;
            final ReadPileup basePileup = alignmentContext.getBasePileup();
            return basePileup.getPileupString(bases != null ? (char) bases[0] : 'N');
        }).saveAsTextFile(outputFile);
    }

    // similar to ReadShard, except there is no underlying ReadSource
    private List<SimpleInterval> getIntervalShards(List<SimpleInterval> intervals, int intervalShardSize) {
        final List<SimpleInterval> shards = new ArrayList<>();
        for (SimpleInterval interval : intervals) {
            int start = interval.getStart();
            while (start <= interval.getEnd()) {
                int end = Math.min(start + intervalShardSize - 1, interval.getEnd());

                final SimpleInterval nextShardInterval = new SimpleInterval(interval.getContig(), start, end);
                shards.add(nextShardInterval);

                start += intervalShardSize;
            }
        }
        return shards;
    }

    @SuppressWarnings("unchecked")
    private static class PileupFunction implements FlatMapFunction2<Iterator<GATKRead>, Iterator<SimpleInterval>, AlignmentContext> {
        private static final long serialVersionUID = 1L;

        private SAMFileHeader header;

        public PileupFunction(SAMFileHeader header) {
            this.header = header;
        }

        @Override
        public Iterable<AlignmentContext> call(Iterator<GATKRead> readIterator,
                                                         Iterator<SimpleInterval> intervalIterator) {
            // get the samples from the read groups
            final Set<String> samples = header.getReadGroups().stream()
                    .map(SAMReadGroupRecord::getSample)
                    .collect(Collectors.toSet());
            LocusIteratorByState libs = new LocusIteratorByState(readIterator, DownsamplingMethod.NONE, true, false, false, samples, header);
            return new IntervalOverlappingIterator<>(libs, ImmutableList.copyOf(intervalIterator), header.getSequenceDictionary());
        }
    }

    private static class ShufflePileupFunction implements FlatMapFunction<Tuple2<SimpleInterval, Iterable<GATKRead>>, AlignmentContext> {
        private static final long serialVersionUID = 1L;

        private SAMFileHeader header;

        public ShufflePileupFunction(SAMFileHeader header) {
            this.header = header;
        }

        @Override
        public Iterable<AlignmentContext> call(Tuple2<SimpleInterval, Iterable<GATKRead>> t) throws Exception {
            SimpleInterval interval = t._1();
            Iterator<GATKRead> readIterator = t._2().iterator();
            // get the samples from the read groups
            final Set<String> samples = header.getReadGroups().stream()
                    .map(SAMReadGroupRecord::getSample)
                    .collect(Collectors.toSet());
            LocusIteratorByState libs = new LocusIteratorByState(readIterator, DownsamplingMethod.NONE, true, false, false, samples, header);
            return new IntervalOverlappingIterator<>(libs, ImmutableList.of(interval), header.getSequenceDictionary());
        }
    }

}
