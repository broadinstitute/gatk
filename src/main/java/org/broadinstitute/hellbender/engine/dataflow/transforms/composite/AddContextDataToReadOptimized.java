package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.DoFnWLog;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipListOneContig;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;


/**
 * A sequence of steps that progressively fill in the ContextShard.
 *
 * Initially we have only the interval and variants.
 * subdivideAndFillReads then loads the reads (and splits into smaller shards)
 * fillContext then adds the reference bases.
 *
 * In the end we have reads, variants, and reference bases.
 *
 * Or if you prefer less typing, you can use add.
 */
public class AddContextDataToReadOptimized {

    /**
     * Immutable holder class.
     * This class enforces no consistency guarantees. The intent of course
     * is for the reads and variants to be in the interval, but that's up the callers.
     */
    public static class ContextShard implements Serializable {
        private static final long serialVersionUID = 1L;

        // the interval covered by this shard
        public final SimpleInterval interval;
        // variants that overlap with the shard.
        public final IntervalsSkipListOneContig<Variant> variants;
        // reads that start in the shard
        public final ArrayList<GATKRead> reads;
        // variants and reference for the particular read at the same index as this element.
        public final ArrayList<ReadContextData> readContext;

        public ContextShard(SimpleInterval interval) {
            this.interval = interval;
            this.variants = null;
            this.reads = new ArrayList<>();
            this.readContext = new ArrayList<>();
        }

        /**
         * Careful: this ctor takes ownership of the passed reads and ReadContextData array.
         * Do not modify them after this call (ideally don't even keep a reference to them).
         */
        private ContextShard(SimpleInterval interval, IntervalsSkipListOneContig<Variant> variants, final ArrayList<GATKRead> reads, final ArrayList<ReadContextData> readContext) {
            this.interval = interval;
            this.variants = variants;
            this.reads = reads;
            this.readContext = readContext;
        }

        /**
         * create a new shard, keeping only the variants that overlap
         * with the new interval. Reads etc are unchanged.
         */
        public ContextShard split(SimpleInterval newInterval) {
            final IntervalsSkipListOneContig<Variant> newVariants;
            if (null==variants) {
                newVariants = null;
            } else {
                newVariants = new IntervalsSkipListOneContig<>( variants.getOverlapping(newInterval) );
            }
            return new ContextShard(newInterval, newVariants, reads, readContext);
        }

        /**
         * creates a new shard, adding the specified variants.
         * Note that readContext is unchanged (including the variants it may refer to).
         */
        public ContextShard withVariants(ArrayList<Variant> newVariants) {
            return new ContextShard(this.interval, new IntervalsSkipListOneContig<>(newVariants), reads, readContext);
        }

        /**
         * creates a new shard, adding the specified reads.
         * Careful: this call takes ownership of the passed read array. So you can't modify it after this call.
         */
        public ContextShard withReads(ArrayList<GATKRead> newReads) {
            return new ContextShard(this.interval, this.variants, newReads, readContext);
        }

        /**
         * creates a new shard, adding the specified read context and *removing the variants*.
         * Careful: this call takes ownership of the passed ReadContextData array. So you can't modify it after this call.
         */
        public ContextShard withReadContext(ArrayList<ReadContextData> newReadContext) {
            return new ContextShard(interval, null, reads, newReadContext);
        }

        /**
         * Returns the variants that overlap the query interval, in start-position order.
         */
        public ArrayList<Variant> variantsOverlapping(SimpleInterval interval) {
            return variants.getOverlapping(interval);
        }

    }

    /**
     * Takes the variants, groups them into shards of size "bigShardSize" at the client, then
     * runs subdivideAndFillReads then fillContext. See their documentation for more information.
     */
    public static PCollection<ContextShard> add(
            Pipeline pipeline,
            List<SimpleInterval> intervalsOfInterest, List<Variant> variants,
            String bam, final ReadFilter optFilter,
            final ReferenceDataflowSource refSource) throws IOException {

        return AddContextDataToReadOptimized.add(pipeline, intervalsOfInterest, 1_000_000, variants,
                bam, 5000, 1000, optFilter,
                refSource);
    }

    /**
     * Takes the variants, groups them into shards of size "bigShardSize" at the client, then
     * runs subdivideAndFillReads then fillContext. See their documentation for more information.
     */
    public static PCollection<ContextShard> add(
            Pipeline pipeline,
            List<SimpleInterval> intervalsOfInterest, int bigShardSize, List<Variant> variants,
            String bam, int outputShardSize, int margin, final ReadFilter optFilter,
            final ReferenceDataflowSource refSource
    ) throws IOException {

        List<SimpleInterval> shardedIntervals = IntervalUtils.cutToShards(intervalsOfInterest, bigShardSize);
        ArrayList<AddContextDataToReadOptimized.ContextShard> shards = AddContextDataToReadOptimized.fillVariants(shardedIntervals, variants, margin);
        PCollection<AddContextDataToReadOptimized.ContextShard> shardsPCol = pipeline.apply(Create.of(shards));
        return shardsPCol
                // big shards of variants -> smaller shards with variants, reads. We take the opportunity to filter the reads as close to the source as possible.
                .apply(ParDo.named("subdivideAndFillReads").of(AddContextDataToReadOptimized.subdivideAndFillReads(bam, outputShardSize, margin, optFilter)))
                        // add ref bases to the shards.
                .apply(ParDo.named("fillContext").of(AddContextDataToReadOptimized.fillContext(refSource)));

    }

    /**
     * Fill in reads that start in the given shard, and subshard the output to the requested size.
     * Signals an error if any read sticks out more than "margin" outside of the big shard.
     * Optionally filters the reads to only keep the ones that satisfy the given predicate.
     *
     * Hadoop paths aren't yet supported.
     *
     *   _   _  ____ _______ ______
     *  | \ | |/ __ \__   __|  ____|
     *  |  \| | |  | | | |  | |__
     *  | . ` | |  | | | |  |  __|
     *  | |\  | |__| | | |  | |____
     *  |_| \_|\____/  |_|  |______| (important):
     *
     * The reads are stored without the header.
     *
     * @param bam URL to the reads
     * @param outputShardSize subsharding size
     * @param margin throw an exception if any read sticks out beyond the margin
     * @param optFilter if specified, only reads that satisfy this will be included.
     * @return A DoFn that acts as described above.
     * @throws IOException
     */
    public static DoFn<ContextShard,ContextShard> subdivideAndFillReads(String bam, int outputShardSize, int margin, final ReadFilter optFilter) throws IOException {
        return new DoFnWLog<ContextShard, ContextShard>("subdivideAndFillReads") {
            private static final long serialVersionUID = 1L;
            private SamReader reader = null;
            private Storage.Objects storageClient = null;

            @Override
            public void startBundle(Context c) throws Exception {
                super.startBundle(c);
                if (BucketUtils.isCloudStorageUrl(bam)) {
                    storageClient = DataflowCommandLineProgram.HellbenderDataflowOptions.Methods.createStorageClient(c.getPipelineOptions().as(DataflowCommandLineProgram.HellbenderDataflowOptions.class));
                    reader = BAMIO.openBAM(storageClient, bam, ValidationStringency.SILENT);
                } else if (BucketUtils.isHadoopUrl(bam)) {
                    throw new RuntimeException("Sorry, Hadoop paths aren't yet supported");
                } else {
                    // read from local file (this only makes sense for the direct runner)
                    reader = SamReaderFactory.make().validationStringency(ValidationStringency.SILENT).open(new File(bam));
                }
            }

            @Override
            public void processElement(ProcessContext c) throws Exception {
                ContextShard shard = c.element();
                ArrayList<SimpleInterval> ints =new ArrayList<>();
                ints.add(shard.interval);
                List<SimpleInterval> subshards = IntervalUtils.cutToShards(ints, outputShardSize);
                int currentSubShardIndex = 0;
                final int lastValidPos = shard.interval.getEnd() + margin;
                final int firstValidPos = shard.interval.getStart() - margin;
                SimpleInterval currentSubShard = subshards.get(currentSubShardIndex);

                try (SAMRecordIterator query = reader.queryOverlapping(shard.interval.getContig(), shard.interval.getStart(), shard.interval.getEnd())) {
                    ArrayList<GATKRead> readsSoFar = new ArrayList<>();

                    while (query.hasNext()) {
                        SAMRecord r = query.next();
                        SAMRecordToGATKReadAdapter g = new SAMRecordToGATKReadAdapter(r);
                        // yes, it'd be a tad faster to check before the wrapping.
                        // But this keeps the code a tad simpler.
                        if (!accept(g, shard.interval)) continue;
                        if (null!=optFilter) {
                            // skip reads that don't pass the filter
                            if (!optFilter.test(g)) continue;
                        }
                        if (!g.isUnmapped()) {
                            // error out if we accept a read that sticks out too far
                            // (margin was too tight, the shard may end up missing relevant variants)
                            if (r.getAlignmentEnd()>lastValidPos) {
                                throw new GATKException("Margin was too tight, a read sticks out by "+(r.getAlignmentEnd()-shard.interval.getEnd())+", going all the way to "+r.getAlignmentEnd());
                            }
                            if (r.getAlignmentStart()<firstValidPos) {
                                throw new GATKException("Margin was too tight, a read starts early by "+(shard.interval.getStart()-r.getAlignmentStart())+", starting at "+r.getAlignmentStart());
                            }
                        }

                        while (currentSubShard.getEnd() < r.getStart()) {
                            if (!readsSoFar.isEmpty()) {
                                // ship this one.
                                ContextShard ret = shard.split(currentSubShard).withReads(readsSoFar);
                                c.output(ret);
                                readsSoFar = new ArrayList<>();
                            }
                            // move to the next shard
                            currentSubShard = subshards.get(++currentSubShardIndex);
                        }
                        // the header slows serialization too much
                        g.setHeader(null);
                        readsSoFar.add(g);
                    }
                    // done reading, ship what we have
                    if (!readsSoFar.isEmpty()) {
                        ContextShard ret = shard.split(currentSubShard).withReads(readsSoFar);
                        c.output(ret);
                    }
                }
            }

            @Override
            public void finishBundle(Context c) throws Exception {
                if (null!=reader) {
                    reader.close();
                }
                super.finishBundle(c);
            }

            private boolean accept(GATKRead r, SimpleInterval region) {
                if (r.isUnmapped()) {
                    return (!r.mateIsUnmapped() && r.getMateStart() >= region.getStart());
                } else {
                    // mapped read.
                    // the query gives us reads that overlap the region. We want only reads that *start* in those intervals.
                    return r.getStart() >= region.getStart();
                }
            }

        };
    };

    /**
     * Given a shard that has reads and variants, query Google Genomics' Reference server and get reference info
     * (including an extra margin on either side), and fill that and the correct variants into readContext.
     */
    public static DoFn<ContextShard, ContextShard> fillContext(final ReferenceDataflowSource refSource) throws IOException {
        return new DoFnWLog<ContextShard, ContextShard>("fillContext") {
            private static final long serialVersionUID = 1L;

            @Override
            public void processElement(ProcessContext c) throws Exception {
                ContextShard shard = c.element();

                // use the function to make sure we get the exact correct amount of reference bases
                int start = Integer.MAX_VALUE;
                int end = Integer.MIN_VALUE;
                SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction = refSource.getReferenceWindowFunction();
                for (GATKRead r : shard.reads) {
                    SimpleInterval readRefs = referenceWindowFunction.apply(r);
                    start = Math.min(readRefs.getStart(), start);
                    end = Math.max(readRefs.getEnd(), end);
                }

                if (start==Integer.MAX_VALUE) {
                    // there are no reads in this shard, so we're going to remove it
                    return;
                }

                SimpleInterval refInterval = new SimpleInterval(shard.interval.getContig(), start, end);
                ReferenceBases refBases = refSource.getReferenceBases(c.getPipelineOptions(), refInterval);

                ArrayList<ReadContextData> readContext = new ArrayList<>();
                for (GATKRead r : shard.reads) {
                    SimpleInterval readInterval = new SimpleInterval(r);
                    ArrayList<Variant> variantsOverlappingThisRead = shard.variantsOverlapping(readInterval);
                    // we pass all the bases. That's better because this way it's just a shared
                    // pointer instead of being an array copy. Downstream processing is fine with having
                    // extra bases (it expects a few, actually).
                    readContext.add(new ReadContextData(refBases, variantsOverlappingThisRead));
                }
                ContextShard ret = shard.withReadContext(readContext);
                c.output(ret);
            }
        };
    }

    /**
     * Goes from full shards to separated, individual (read,context pairs).
     * The advantage is you don't need to know about shards anymore.
     * This comes at the expense of copying the reference data instead of being able to have a pointer to in-shard data.
     */
    public static DoFn<ContextShard,KV<GATKRead,ReadContextData>> flattenShards(SerializableFunction<GATKRead, SimpleInterval> contextFn) {
        return new DoFnWLog<ContextShard, KV<GATKRead, ReadContextData>>("flattenShards") {
            private static final long serialVersionUID = 1L;
            @Override
            public void processElement(ProcessContext c) throws Exception {
                ContextShard shard = c.element();
                for (int i=0; i<shard.reads.size(); i++) {
                    GATKRead read = shard.reads.get(i);
                    ReadContextData rc = shard.readContext.get(i);
                    ReadContextData rcd = new ReadContextData( rc.getOverlappingReferenceBases().getSubset(contextFn.apply(read)), rc.getOverlappingVariants());
                    c.output(KV.of(read,rcd));
                }
            }
        };
    }


    /**
     * Given a list of shards and a list of variants,
     * add each variant to every (shard+margin) that it overlaps.
     *
     * This happens immediately, at the caller.
     */
    public static ArrayList<AddContextDataToReadOptimized.ContextShard> fillVariants(List<SimpleInterval> shardedIntervals, List<Variant> variants, int margin) {
        IntervalsSkipList<Variant> intervals = new IntervalsSkipList<>(variants);
        ArrayList<AddContextDataToReadOptimized.ContextShard> ret = new ArrayList<>();
        for (SimpleInterval s : shardedIntervals) {
            int start = Math.max(s.getStart() - margin, 1);
            int end = s.getEnd() + margin;
            // here it's OK if end is past the contig's boundary, there just won't be any variant there.
            SimpleInterval interval = new SimpleInterval(s.getContig(), start, end);
            ret.add(new AddContextDataToReadOptimized.ContextShard(s).withVariants(intervals.getOverlapping(interval)));
        }
        return ret;
    }

}