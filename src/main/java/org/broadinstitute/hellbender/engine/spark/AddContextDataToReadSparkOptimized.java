package org.broadinstitute.hellbender.engine.spark;

import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.ContextShard;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;


import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;


public final class AddContextDataToReadSparkOptimized implements Serializable {
    private static final long serialVersionUID = 1L;

    // the granularity at which we'll want to assign work and read inputs.
    public static final int bigShardSize = 1_000_000;
    // the granularity at which we want to batch reads,variants,and reference for processing.
    public static final int outputShardSize = 5_000;
    // we make the assumption that all reads that start in a shard
    // will end within "margin" of the shard. This allows us to split the computation
    // across machines.
    public static final int margin = 1000;


    /**
     * Create shards with reads, variants, and reference bases, using default values for shard sizes and margin.
     * See the other methods here for an explanation of the various arguments.
     */
     public static JavaRDD<ContextShard> add(JavaSparkContext ctx, final List<SimpleInterval> intervals,
                                             String bam, final List<GATKVariant> variants, AuthHolder auth,
                                             final ReadFilter optFilter, final ReferenceMultiSource rds) {
        // prepare shards for the intervals of interest
        List<SimpleInterval> shardedIntervals = IntervalUtils.cutToShards(intervals, bigShardSize);
        // add variants
        ArrayList<ContextShard> localShards = AddContextDataToReadSparkOptimized.fillVariants(shardedIntervals, variants, margin);
        // ship to cluster
        JavaRDD<ContextShard> shards = ctx.parallelize(localShards);
        // subdivide, and add reads
        JavaRDD<ContextShard> reads;
        try {
            reads = shards.flatMap(AddContextDataToReadSparkOptimized.subdivideAndFillReads(bam, auth, outputShardSize, margin, optFilter));
        } catch (IOException x) {
            throw new UserException.CouldNotReadInputFile("Couldn't read "+bam+": "+x.getMessage(), x);
        }
        // add reference bases
        reads = reads.map(s -> AddContextDataToReadSparkOptimized.fillContext(rds, s));

        return reads;
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
     * @param bam URL to the reads. CRAM files aren't supported yet (nor are Hadoop files).
     * @param outputShardSize subsharding size
     * @param margin throw an exception if any read sticks out beyond the margin
     * @param optFilter if specified, only reads that satisfy this will be included.
     * @return A FlatMapFunction that acts as described above.
     * @throws IOException
     */
    public static FlatMapFunction<ContextShard,ContextShard> subdivideAndFillReads(String bam, AuthHolder auth, int outputShardSize, int margin, final ReadFilter optFilter) throws IOException {
            return new FlatMapFunction<ContextShard, ContextShard>() {
                private static final long serialVersionUID = 1L;
                @Override
                public Iterable<ContextShard> call(ContextShard contextShard) throws Exception {
                    return new Iterable<ContextShard>() {
                        @Override
                        public Iterator<ContextShard> iterator() {
                            try {
                                return new SubdivideAndFillReadsIterator(bam, auth, outputShardSize, margin, optFilter, contextShard);
                            } catch (Exception x) {
                                throw new RuntimeException(x);
                            }
                        }
                    };
                }
            };
    }

    /**
     * Given a shard that has reads and variants, query Google Genomics' Reference server and get reference info
     * (including an extra margin on either side), and fill that and the correct variants into readContext.
     */
    public static ContextShard fillContext(ReferenceMultiSource refSource, ContextShard shard) {
        if (null==shard) return null;
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
            return null;
        }

        SimpleInterval refInterval = new SimpleInterval(shard.interval.getContig(), start, end);
        ReferenceBases refBases;
        try {
            refBases = refSource.getReferenceBases(null, refInterval);
        } catch (IOException x) {
            throw new GATKException("Unable to read the reference");
        }

        ArrayList<ReadContextData> readContext = new ArrayList<>();
        for (GATKRead r : shard.reads) {
            SimpleInterval readInterval = new SimpleInterval(r);
            List<GATKVariant> variantsOverlappingThisRead = shard.variantsOverlapping(readInterval);
            // we pass all the bases. That's better because this way it's just a shared
            // pointer instead of being an array copy. Downstream processing is fine with having
            // extra bases (it expects a few, actually).
            readContext.add(new ReadContextData(refBases, variantsOverlappingThisRead));
        }
        return shard.withReadContext(readContext);
    }

    protected static final class SubdivideAndFillReadsIterator implements Iterator<ContextShard>, Serializable {
        private static final long serialVersionUID = 1L;
        private static final Logger log = LogManager.getLogger(SubdivideAndFillReadsIterator.class);
        private static final int maxReadsPerShard = 10_000;
        private final transient SamReader reader;
        private final transient SAMRecordIterator query;
        private final ContextShard shard;
        private final String bam;
        private final ReadFilter optFilter;
        private final int lastValidPos;
        private final int firstValidPos;
        private final List<SimpleInterval> subshards;
        private int currentSubShardIndex;
        private ArrayList<GATKRead> readsSoFar = new ArrayList<>();
        private Storage.Objects storageClient = null;
        private SimpleInterval currentSubShard;
        private ContextShard nextOutput = null;
        private boolean readerClosed = false;

        public SubdivideAndFillReadsIterator(String bam, AuthHolder auth, int outputShardSize, int margin, final ReadFilter optFilter, ContextShard shard) throws IOException, GeneralSecurityException, ClassNotFoundException {
            this.bam = bam;
            this.shard = shard;
            this.optFilter = optFilter;
            // it's OK if this goes beyond the contig boundaries.
            lastValidPos = shard.interval.getEnd() + margin;
            firstValidPos = Math.max(shard.interval.getStart() - margin, 1);
            ArrayList<SimpleInterval> ints =new ArrayList<>();
            ints.add(shard.interval);
            subshards = IntervalUtils.cutToShards(ints, outputShardSize);
            currentSubShardIndex = 0;
            currentSubShard = subshards.get(currentSubShardIndex);

            if (BucketUtils.isCloudStorageUrl(bam)) {
                storageClient = auth.makeStorageClient();
                reader = BAMIO.openBAM(storageClient, bam, ValidationStringency.SILENT);
            } else if (BucketUtils.isHadoopUrl(bam)) {
                throw new RuntimeException("Sorry, Hadoop paths aren't yet supported");
            } else {
                // read from local file (this only makes sense if every worker sees the same thing, e.g. if we're running locally)
                reader = SamReaderFactory.make().validationStringency(ValidationStringency.SILENT).open(new File(bam));
            }
            query = reader.queryOverlapping(shard.interval.getContig(), shard.interval.getStart(), shard.interval.getEnd());

        }

        @Override
        public boolean hasNext() {
            if (null==nextOutput) nextOutput = tryNext();
            return (null!=nextOutput);
        }

        @Override
        public ContextShard next() {
            if (null==nextOutput) nextOutput = tryNext();
            if (null==nextOutput) throw new NoSuchElementException();
            ContextShard ret =  nextOutput;
            nextOutput = null;
            return ret;
        }

        // returns the next shard, if any. Otherwise, closes the reader and returns null.
        private ContextShard tryNext() {
            if (readerClosed) return null;
            ContextShard ret = null;
            while (query.hasNext()) {

                SAMRecord r = query.next();
                SAMRecordToGATKReadAdapter g = new SAMRecordToGATKReadAdapter(r);
                // yes, it'd be a tad faster to check before the wrapping.
                // But this keeps the code a tad simpler.
                if (!accept(g, shard.interval)) continue;
                // we expect it's in, but what if we're wrong? Better to detect it than return a wrong answer.
                throwIfOutsideMargin(g, r);

                while (currentSubShard.getEnd() < r.getStart()) {
                    if (!readsSoFar.isEmpty()) {
                        // ship this one.
                        ret = shard.split(currentSubShard).withReads(readsSoFar);
                        readsSoFar = new ArrayList<>();
                    }
                    // move to the next shard.
                    // This cannot overflow because there is no read that start past the final shard ending,
                    // so the while condition above will be false on the last shard.
                    currentSubShard = subshards.get(++currentSubShardIndex);
                }
                // the header slows serialization too much
                g.setHeader(null);
                readsSoFar.add(g);
                if (readsSoFar.size()>=maxReadsPerShard) {
                    // ship this one.
                    log.info("Too many reads in this shard, splitting it."+readsSoFar.size());
                    int currentSubShardEnd;
                    if (g.isUnmapped()) {
                        if (!g.mateIsUnmapped()) {
                            currentSubShardEnd = g.getMateStart() + margin;
                        } else {
                            throw new GATKException.ShouldNeverReachHereException("How did an unmapped read make it to here? "+g.toString());
                        }
                    } else {
                        currentSubShardEnd = g.getStart() + margin;
                    }
                    // we grow the interval by "margin" to make sure we get all the variants we need.
                    // Since we already assume that margin has that property, we're good to go.
                    SimpleInterval thisInterval = new SimpleInterval(currentSubShard.getContig(), currentSubShard.getStart(), currentSubShardEnd);
                    ret = shard.split(thisInterval).withReads(readsSoFar);
                    readsSoFar = new ArrayList<>();
                    // do not advance currentSubShard, we're still technically in the same one.
                }

                if (null!=ret) {
                    // we have something to ship, and we're done handing the read we just got out of the query iterator.
                    return ret;
                }
            }
            // done reading, ship what we have
            readerClosed = true;
            query.close();
            try {
                if (null!=reader) {
                    reader.close();
                }
            } catch (IOException x) {
                throw new GATKException.ShouldNeverReachHereException("IOException when closing the BAM file reader for "+bam);
            }
            if (!readsSoFar.isEmpty()) {
                return shard.split(currentSubShard).withReads(readsSoFar);
            }
            return null;
        }

        // check mapping and optFilter
        private boolean accept(GATKRead r, SimpleInterval region) {
            boolean ret;
            if (r.isUnmapped()) {
                ret = (!r.mateIsUnmapped() && r.getMateStart() >= region.getStart());
            } else {
                // mapped read.
                // the query gives us reads that overlap the region. We want only reads that *start* in those intervals.
                ret = r.getStart() >= region.getStart();
            }
            // apply optional filter, if specified
            return ret && (null==optFilter || optFilter.test(r));
        }


        private void throwIfOutsideMargin(SAMRecordToGATKReadAdapter g, SAMRecord r) {
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
        }

    }



    /**
     * Given a list of shards and a list of variants,
     * add each variant to every (shard+margin) that it overlaps.
     *
     * This happens immediately, at the caller.
     */
    public static ArrayList<ContextShard> fillVariants(List<SimpleInterval> shardedIntervals, List<GATKVariant> variants, int margin) {
        IntervalsSkipList<GATKVariant> intervals = new IntervalsSkipList<>(variants);
        ArrayList<ContextShard> ret = new ArrayList<>();
        for (SimpleInterval s : shardedIntervals) {
            int start = Math.max(s.getStart() - margin, 1);
            int end = s.getEnd() + margin;
            // here it's OK if end is past the contig's boundary, there just won't be any variant there.
            SimpleInterval expandedInterval = new SimpleInterval(s.getContig(), start, end);
            // the next ContextShard has interval s because we want it to contain all reads that start in s.
            // We give it all variants that overlap the expanded interval in order to make sure we include
            // all the variants that overlap with the reads of interest.
            //
            // Graphically:
            // |------- s --------|
            //--------expandedInterval------------------|
            //            |-- a read starting in s --|
            //                           |--- a variant overlapping the read ---|
            //
            // Since the read's length is less than margin, we know that by including all the variants that overlap
            // with the expanded interval we are also including all the variants that overlap with all the reads in this shard.
            ret.add(new ContextShard(s).withVariants(intervals.getOverlapping(expandedInterval)));
        }
        return ret;
    }

}