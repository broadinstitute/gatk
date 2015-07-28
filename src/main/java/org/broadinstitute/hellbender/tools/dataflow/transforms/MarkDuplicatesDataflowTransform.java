package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.cloud.dataflow.sdk.coders.CoderException;
import com.google.cloud.dataflow.sdk.coders.CustomCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.coders.StringUtf8Coder;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.Flatten;
import com.google.cloud.dataflow.sdk.transforms.GroupByKey;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.Partition;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionList;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimaps;
import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.tools.dataflow.MarkDuplicatesDataflow;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * The main trasform for MarkDuplicates. Takes reads, returned reads with some reads marked as duplicates.
 */
public final class MarkDuplicatesDataflowTransform extends PTransform<PCollection<GATKRead>, PCollection<GATKRead>> {
    private static final long serialVersionUID = 1l;

    //private vs non-primary alignments
    private enum ReadsPartition {
        PRIMARY, NOT_PRIMARY
    }

    //Bases below this quality will not be included in picking the best read from a set of duplicates.
    private static final int MIN_BASE_QUAL = 15;



    /**
     * Struct-like class to store information about the paired reads for mark duplicates.
     */
    private static final class PairedEnds {
        private GATKRead first, second;

        PairedEnds(final GATKRead first) {
            this.first = first;
        }

        public static PairedEnds of(final GATKRead first) {
            return new PairedEnds(first);
        }

        public PairedEnds and(final GATKRead second) {
            if (second != null &&
                ReadUtils.getStrandedUnclippedStart(first) > ReadUtils.getStrandedUnclippedStart(second)) {

                this.second = this.first;
                this.first = second;
            } else {
                this.second = second;
            }
            return this;
        }

        public String key(final SAMFileHeader header) {
            return MarkDuplicatesReadsKey.keyForPairedEnds(header, first, second);
        }

        public GATKRead first() {
            return first;
        }

        public GATKRead second() {
            return second;
        }

        public int score() {
            return MarkDuplicatesDataflowTransform.score(first) + MarkDuplicatesDataflowTransform.score(second);
        }
    }

    /**
     * Special coder for the PairedEnds class.
     */
    private static final class PairedEndsCoder extends CustomCoder<PairedEnds> {
        private static final long serialVersionUID = 1L;

        private static final CustomCoder<GATKRead> readCoder = new GATKReadCoder();

        @Override
        public void encode( PairedEnds value, OutputStream outStream, Context context ) throws CoderException, IOException {
            if ( value == null || value.first() == null ) {
                throw new IOException("nothing to encode");
            }
            final boolean isCompletePair = value.second() != null ;
            SerializableCoder.of(Boolean.class).encode(isCompletePair, outStream, context);

            readCoder.encode(value.first, outStream, context);
            if ( isCompletePair ) {
                readCoder.encode(value.second, outStream, context);
            }
        }

        @Override
        public PairedEnds decode( InputStream inStream, Context context ) throws CoderException, IOException {
            final boolean isCompletePair = SerializableCoder.of(Boolean.class).decode(inStream, context);

            final PairedEnds pairedEnds = PairedEnds.of(readCoder.decode(inStream, context));
            if ( isCompletePair ) {
                pairedEnds.and(readCoder.decode(inStream, context));
            }
            return pairedEnds;
        }
    }

    // TODO: extract this into an independent class, and unify with other comparators in the codebase
    private static final class CoordinateOrder implements Comparator<GATKRead>, Serializable {
        private static final long serialVersionUID = 1l;
        private final SAMFileHeader header;

        public CoordinateOrder(final SAMFileHeader header) {
            this.header = header;
        }

        @Override
        public int compare(final GATKRead lhs, final GATKRead rhs) {
            if (rhs == lhs) return 0; //shortcut

            final int res1 = Integer.compare(ReadUtils.getReferenceIndex(lhs, header), ReadUtils.getReferenceIndex(rhs, header));
            if (res1 != 0) return res1;

            final int res2 = Long.compare(lhs.getStart(), rhs.getStart());
            if (res2 != 0) return res2;

            final int res3 = Boolean.compare(lhs.isDuplicate(), rhs.isDuplicate());
            if (res3 != 0) return res3;

            final int res4 = Boolean.compare(lhs.failsVendorQualityCheck(), rhs.failsVendorQualityCheck());
            if (res4 != 0) return res4;

            final int res5 = Boolean.compare(lhs.isPaired(), rhs.isPaired());
            if (res5 != 0) return res5;

            final int res6 = Boolean.compare(lhs.isProperlyPaired(), rhs.isProperlyPaired());
            if (res6 != 0) return res6;

            final int res7 = Boolean.compare(lhs.isFirstOfPair(), rhs.isFirstOfPair());
            if (res7 != 0) return res7;

            final int res8 = Boolean.compare(lhs.isSecondaryAlignment(), rhs.isSecondaryAlignment());
            if (res8 != 0) return res8;

            final int res9 = Boolean.compare(lhs.isSupplementaryAlignment(), rhs.isSupplementaryAlignment());
            if (res9 != 0) return res9;

            final int res10 = Integer.compare(lhs.getMappingQuality(), rhs.getMappingQuality());
            if (res10 != 0) return res10;

            final int res11 = Integer.compare(ReadUtils.getMateReferenceIndex(lhs, header), ReadUtils.getMateReferenceIndex(rhs, header));
            if (res11 != 0) return res11;

            final int res12 = Long.compare(lhs.getMateStart(), rhs.getMateStart());
            return res12 ;
        }
    }

    /**
     *  The header for all reads processed.
     */
    private final PCollectionView<SAMFileHeader> header;

    public MarkDuplicatesDataflowTransform( final PCollectionView<SAMFileHeader> header ) {
        this.header = header;
    }

    @Override
    public PCollection<GATKRead> apply( final PCollection<GATKRead> preads ) {
        final PCollectionList<GATKRead> readsPartitioned = partitionReadsByPrimaryNonPrimaryAlignment(preads);

        final PCollection<GATKRead> fragments = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
        final PCollection<GATKRead> fragmentsTransformed = transformFragments(header, fragments);

        final PCollection<GATKRead> pairs = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
        final PCollection<GATKRead> pairsTransformed = transformReads(header, pairs);

        //no work on those
        final PCollection<GATKRead> not_primary = readsPartitioned.get(ReadsPartition.NOT_PRIMARY.ordinal());

        return PCollectionList.of(fragmentsTransformed).and(pairsTransformed).and(not_primary).apply(Flatten.<GATKRead>pCollections());
    }

    private static PCollectionList<GATKRead> partitionReadsByPrimaryNonPrimaryAlignment(PCollection<GATKRead> preads) {
        return preads.apply(Partition.of(2, new Partition.PartitionFn<GATKRead>() {
            private static final long serialVersionUID = 1l;

            @Override
            public int partitionFor(final GATKRead read, final int n) {
                return read.isSecondaryAlignment() || read.isSupplementaryAlignment() || read.isUnmapped() ? ReadsPartition.NOT_PRIMARY.ordinal() : ReadsPartition.PRIMARY.ordinal();
            }
        }));
    }
    /**
     * (1) Reads are grouped by read group and read name.
     * (2) Reads are then paired together as follows:
     *   (a) The remaining reads (one per fragment key) are coordinate-sorted and paired
     *       consecutively together and emitted.
     *   (b) If a read is leftover, it is emitted, unmodified, as an unpaired end.
     * (3) Paired ends are grouped by a similar key as above but using both reads.
     * (4) Any unpaired end is emitted, unmodified.
     * (5) The remained paired ends are scored and all but the highest scoring are marked as
     *     duplicates. Both reads in the pair are emitted.
     */
    private PCollection<GATKRead> transformReads(final PCollectionView<SAMFileHeader> headerPcolView, final PCollection<GATKRead> pairs) {
        final PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> makeKeysForPairs = makeKeysForPairs(headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<KV<String, PairedEnds>>> markGroupedDuplicatePairs = markGroupedDuplicatePairs(headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<GATKRead>> markPairedEnds = markPairedEnds();

        return pairs
            .apply(makeKeysForPairs)
            .apply(GroupByKey.<String, GATKRead>create())
            .apply(markGroupedDuplicatePairs)
            .setCoder(KvCoder.of(StringUtf8Coder.of(), new PairedEndsCoder()))
            .apply(GroupByKey.<String, PairedEnds>create())
            .apply(markPairedEnds);
    }

    /**
     * Makes keys for read pairs. To be grouped by in the next step.
     */
    private PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> makeKeysForPairs(final PCollectionView<SAMFileHeader> headerPcolView) {
        return ParDo
            .named("make keys for pairs")
            .withSideInputs(headerPcolView)
            .of(new DoFn<GATKRead, KV<String, GATKRead>>() {
                private static final long serialVersionUID = 1l;

                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    final GATKRead record = context.element();
                    if (ReadUtils.readHasMappedMate(record)) {
                        final SAMFileHeader h = context.sideInput(headerPcolView);
                        final String key = MarkDuplicatesReadsKey.keyForPair(h, record);
                        final KV<String, GATKRead> kv = KV.of(key, record);
                        context.output(kv);
                    }
                }
            });
    }

    private PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<KV<String, PairedEnds>>> markGroupedDuplicatePairs(final PCollectionView<SAMFileHeader> headerPcolView) {
        return ParDo
            .named("pair ends")
            .withSideInputs(headerPcolView)
            .of(new DoFn<KV<String, Iterable<GATKRead>>, KV<String, PairedEnds>>() {
                private static final long serialVersionUID = 1L;
                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    final SAMFileHeader header = context.sideInput(headerPcolView);
                    final List<GATKRead> sorted = Lists.newArrayList(context.element().getValue());
                    sorted.sort(new CoordinateOrder(header));
                    PairedEnds pair = null;
                    //Records are sorted, we iterate over them and pair them up.
                    for (final GATKRead record : sorted) {
                        if (pair == null) {                                //first in pair
                            pair = PairedEnds.of(record);
                        } else {                                           //second in pair
                            pair.and(record);
                            context.output(KV.of(pair.key(header), pair));
                            pair = null;                                   //back to first
                        }
                    }
                    if (pair != null) {                                    //left over read
                        context.output(KV.of(pair.key(header), pair));
                    }
                }
            });
    }


    private PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<GATKRead>> markPairedEnds() {
        return ParDo
            .named("mark paired ends")
            .of(new DoFn<KV<String, Iterable<PairedEnds>>, GATKRead>() {
                private static final long serialVersionUID = 1l;

                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    final ImmutableListMultimap<Boolean, PairedEnds> paired = Multimaps.index(context.element().getValue(), pair -> pair.second() != null);

                    // As in Picard, unpaired ends left alone.
                    for (final PairedEnds pair : paired.get(false)) {
                        context.output(pair.first());
                    }

                    //order by score
                    final List<PairedEnds> scored = Ordering.natural().reverse().onResultOf((PairedEnds pair) -> pair.score()).immutableSortedCopy(paired.get(true));
                    final PairedEnds best = Iterables.getFirst(scored, null);
                    if (best != null) {
                        context.output(best.first());
                        context.output(best.second());
                    }
                    //Mark everyone who's not best as a duplicate
                    for (final PairedEnds pair : Iterables.skip(scored, 1)) {
                        GATKRead record = pair.first();
                        record.setIsDuplicate(true);
                        context.output(record);

                        record = pair.second();
                        record.setIsDuplicate(true);
                        context.output(record);
                    }
                }
            });
    }


    /**
     * Takes the reads,
     * group them by library, contig, position and orientation,
     * within each group
     *   (a) if there are only fragments, mark all but the highest scoring as duplicates, or,
     *   (b) if at least one is marked as paired, mark all fragments as duplicates.
     *  Note: Emit only the fragments, as the paired reads are handled separately.
     */
    PCollection<GATKRead> transformFragments(final PCollectionView<SAMFileHeader> headerPcolView, final PCollection<GATKRead> fragments) {
        final PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> makeKeysForFragments =  makeKeysForFragments(headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<GATKRead>> markGroupedDuplicateFragments = markGroupedDuplicateFragments();
        return fragments
            .apply(makeKeysForFragments)
            .apply(GroupByKey.<String, GATKRead>create())
            .apply(markGroupedDuplicateFragments);//no need to set up coder for Read (uses GenericJsonCoder)
    }

    private PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<GATKRead>> markGroupedDuplicateFragments() {
        return ParDo.named("mark dups")
            .of(new DoFn<KV<String, Iterable<GATKRead>>, GATKRead>() {
                private static final long serialVersionUID = 1l;

                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    //split reads by paired vs unpaired
                    final Map<Boolean, List<GATKRead>> byPairing = StreamSupport.stream(context.element().getValue().spliterator(), false).collect(Collectors.partitioningBy(
                        read -> ReadUtils.readHasMappedMate(read)
                    ));

                    // Note the we emit only fragments from this mapper.
                    if (byPairing.get(true).isEmpty()) {
                        // There are no paired reads, mark all but the highest scoring fragment as duplicate.
                        final List<GATKRead> frags = Ordering.natural().reverse().onResultOf((GATKRead read) -> score(read)).immutableSortedCopy(byPairing.get(false));
                        if (!frags.isEmpty()) {
                            context.output(frags.get(0));                         //highest score - just emit
                            for (final GATKRead record : Iterables.skip(frags, 1)) {  //lower   scores - mark as dups and emit
                                record.setIsDuplicate(true);
                                context.output(record);
                            }
                        }
                    } else {
                        // There are paired ends so we mark all fragments as duplicates.
                        for (final GATKRead record : byPairing.get(false)) {
                            record.setIsDuplicate(true);
                            context.output(record);
                        }
                    }
                }
            });
    }

    /**
     * How to assign a score to the read in MarkDuplicates (so that we pick the best one to be the non-duplicate).
     */
    //Note: copied from htsjdk.samtools.DuplicateScoringStrategy
    private static int score(final GATKRead record) {
        if (record == null) {
            return 0;
        } else {
            int sum = 0;
            for ( byte b : record.getBaseQualities() ) {
                int i = (int)b;
                if ( i >= MIN_BASE_QUAL ) {
                    sum += i;
                }
            }
            return sum;
        }
    }


    /**
     * Groups reads by keys - keys are tuples of (library, contig, position, orientation).
     */
    PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> makeKeysForFragments(final PCollectionView<SAMFileHeader> headerPcolView) {
        return ParDo
            .named("make keys for reads")
            .withSideInputs(headerPcolView)
            .of(new DoFn<GATKRead, KV<String, GATKRead>>() {
                private static final long serialVersionUID = 1L;
                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    final GATKRead record = context.element();
                    record.setIsDuplicate(false);
                    final SAMFileHeader h = context.sideInput(headerPcolView);
                    final String key = MarkDuplicatesReadsKey.keyForFragment(h, record);
                    final KV<String, GATKRead> kv = KV.of(key, record);
                    context.output(kv);
                }
            });
    }
}
