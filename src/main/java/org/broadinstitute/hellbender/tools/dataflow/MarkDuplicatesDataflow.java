package org.broadinstitute.hellbender.tools.dataflow;

import com.google.api.client.json.GenericJson;
import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.*;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.transforms.Partition.PartitionFn;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionList;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.cloud.genomics.dataflow.coders.GenericJsonCoder;
import com.google.common.collect.*;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.ReadsSource;
import org.broadinstitute.hellbender.tools.dataflow.transforms.MarkDuplicatesReadsKey;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.tools.dataflow.GenomicsReadUtils.getAlignmentStart;
import static org.broadinstitute.hellbender.tools.dataflow.GenomicsReadUtils.getReferenceName;
import static org.broadinstitute.hellbender.tools.dataflow.GenomicsReadUtils.isPrimaryAlignment;

@CommandLineProgramProperties(
        usage="Marks duplicates on dataflow",
        usageShort="Mark Duplicates",
        programGroup = DataFlowProgramGroup.class)
public final class MarkDuplicatesDataflow extends DataflowCommandLineProgram {
    private static final long serialVersionUID = 1L;

    //Bases below this quality will not be included in picking the best read from a set of duplicates.
    private static final int MIN_BASE_QUAL = 15;

    @Argument(doc="a prefix for the dataflow output files", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "uri for the input bam, either a local file path or a gs:// bucket path",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false)
    protected String bam;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Override
    protected void setupPipeline(final Pipeline pipeline) {
        final ReadsSource readsSource = new ReadsSource(bam, pipeline);
        final SAMFileHeader header = readsSource.getHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(sequenceDictionary):
                getAllIntervalsForReference(sequenceDictionary);

        final PCollectionView<SAMFileHeader> headerPcolView = pipeline.apply(Create.of(header)).apply(View.<SAMFileHeader>asSingleton());

        final PCollection<Read> preads = readsSource.getReadPCollection(intervals);
        final PCollection<Read> results = preads.apply(new MarkDuplicatesDataflowTransform(headerPcolView));

        final PCollection<String> pstrings = results.apply(DataflowUtils.convertToString());
        pstrings.apply(TextIO.Write.to(outputFile));
    }

    /**
     * The main trasform for MarkDuplicates. Takes reads, returned reads with some reads marked as duplicates.
     */
    private static final class MarkDuplicatesDataflowTransform extends PTransform<PCollection<Read>, PCollection<Read>> {

        private static final long serialVersionUID = 1l;

        //private vs non-primary alignments
        private enum ReadsPartition {
            PRIMARY, NOT_PRIMARY
        }

        /**
         *  The header for all reads processed.
        */
        private final PCollectionView<SAMFileHeader> header;

        public MarkDuplicatesDataflowTransform( final PCollectionView<SAMFileHeader> header ) {
            this.header = header;
        }

        @Override
        public PCollection<Read> apply( final PCollection<Read> preads ) {
            final PCollectionList<Read> readsPartitioned = partitionReadsByPrimaryNonPrimaryAlignment(preads);

            final PCollection<Read> fragments = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
            final PCollection<Read> fragmentsTransformed = transformFragments(header, fragments);

            final PCollection<Read> pairs = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
            final PCollection<Read> pairsTransformed = transformReads(header, pairs);

            //no work on those
            final PCollection<Read> not_primary = readsPartitioned.get(ReadsPartition.NOT_PRIMARY.ordinal());

            return PCollectionList.of(fragmentsTransformed).and(pairsTransformed).and(not_primary).apply(Flatten.<Read>pCollections());
        }

        private static PCollectionList<Read> partitionReadsByPrimaryNonPrimaryAlignment(PCollection<Read> preads) {
            return preads.apply(Partition.of(2, new PartitionFn<Read>() {
                private static final long serialVersionUID = 1l;

                @Override
                public int partitionFor(final Read read, final int n) {
                    return isPrimaryAlignment(read) ? ReadsPartition.PRIMARY.ordinal() : ReadsPartition.NOT_PRIMARY.ordinal();
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
        private PCollection<Read> transformReads(final PCollectionView<SAMFileHeader> headerPcolView, final PCollection<Read> pairs) {
            final PTransform<PCollection<? extends Read>, PCollection<KV<String, Read>>> makeKeysForPairs = makeKeysForPairs(headerPcolView);
            final PTransform<PCollection<? extends KV<String, Iterable<Read>>>, PCollection<KV<String, PairedEnds>>> markGroupedDuplicatePairs = markGroupedDuplicatePairs(headerPcolView);
            final PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<Read>> markPairedEnds = markPairedEnds();

            return pairs
                    .apply(makeKeysForPairs)
                    .apply(GroupByKey.<String, Read>create())
                    .apply(markGroupedDuplicatePairs)
                    .setCoder(KvCoder.of(StringUtf8Coder.of(), PairedEndsCoder.of()))
                    .apply(GroupByKey.<String, PairedEnds>create())
                    .apply(markPairedEnds);
        }

        /**
         * Makes keys for read pairs. To be grouped by in the next step.
         */
        private PTransform<PCollection<? extends Read>, PCollection<KV<String, Read>>> makeKeysForPairs(final PCollectionView<SAMFileHeader> headerPcolView) {
            return ParDo
                    .named("make keys for pairs")
                    .withSideInputs(headerPcolView)
                    .of(new DoFn<Read, KV<String, Read>>() {
                        private static final long serialVersionUID = 1l;

                        @Override
                        public void processElement(final ProcessContext context) throws Exception {
                            final Read record = context.element();
                            if (GenomicsReadUtils.isPaired(record)) {
                                final SAMFileHeader h = context.sideInput(headerPcolView);
                                final String key = MarkDuplicatesReadsKey.keyForPair(h, record);
                                final KV<String, Read> kv = KV.of(key, record);
                                context.output(kv);
                            }
                        }
                    });
        }

        private PTransform<PCollection<? extends KV<String, Iterable<Read>>>, PCollection<KV<String, PairedEnds>>> markGroupedDuplicatePairs(final PCollectionView<SAMFileHeader> headerPcolView) {
            return ParDo
                    .named("pair ends")
                    .withSideInputs(headerPcolView)
                    .of(new DoFn<KV<String, Iterable<Read>>, KV<String, PairedEnds>>() {
                        private static final long serialVersionUID = 1L;
                        @Override
                        public void processElement(final ProcessContext context) throws Exception {
                            final SAMFileHeader header = context.sideInput(headerPcolView);
                            final List<Read> sorted = Lists.newArrayList(context.element().getValue());
                            sorted.sort(new CoordinateOrder(header));
                            PairedEnds pair = null;
                            //Records are sorted, we iterate over them and pair them up.
                            for (final Read record : sorted) {
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


        private PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<Read>> markPairedEnds() {
            return ParDo
                    .named("mark paired ends")
                    .of(new DoFn<KV<String, Iterable<PairedEnds>>, Read>() {
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
                                Read record = pair.first();
                                record.setDuplicateFragment(true);
                                context.output(record);

                                record = pair.second();
                                record.setDuplicateFragment(true);
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
        PCollection<Read> transformFragments(final PCollectionView<SAMFileHeader> headerPcolView, final PCollection<Read> fragments) {
            final PTransform<PCollection<? extends Read>, PCollection<KV<String, Read>>> makeKeysForFragments =  makeKeysForFragments(headerPcolView);
            final PTransform<PCollection<? extends KV<String, Iterable<Read>>>, PCollection<Read>> markGroupedDuplicateFragments = markGroupedDuplicateFragments();
            return fragments
                    .apply(makeKeysForFragments)
                    .apply(GroupByKey.<String, Read>create())
                    .apply(markGroupedDuplicateFragments);//no need to set up coder for Read (uses GenericJsonCoder)
        }

        private PTransform<PCollection<? extends KV<String, Iterable<Read>>>, PCollection<Read>> markGroupedDuplicateFragments() {
            return ParDo.named("mark dups")
                    .of(new DoFn<KV<String, Iterable<Read>>, Read>() {
                        private static final long serialVersionUID = 1l;

                        @Override
                        public void processElement(final ProcessContext context) throws Exception {
                            //split reads by paired vs unpaired
                            final Map<Boolean, List<Read>> byPairing = StreamSupport.stream(context.element().getValue().spliterator(), false).collect(Collectors.partitioningBy(
                                    read -> GenomicsReadUtils.isPaired(read)
                            ));

                            // Note the we emit only fragments from this mapper.
                            if (byPairing.get(true).isEmpty()) {
                                // There are no paired reads, mark all but the highest scoring fragment as duplicate.
                                final List<Read> frags = Ordering.natural().reverse().onResultOf((Read read) -> score(read)).immutableSortedCopy(byPairing.get(false));
                                if (!frags.isEmpty()) {
                                    context.output(frags.get(0));                         //highest score - just emit
                                    for (final Read record : Iterables.skip(frags, 1)) {  //lower   scores - mark as dups and emit
                                        record.setDuplicateFragment(true);
                                        context.output(record);
                                    }
                                }
                            } else {
                                // There are paired ends so we mark all fragments as duplicates.
                                for (final Read record : byPairing.get(false)) {
                                    record.setDuplicateFragment(true);
                                    context.output(record);
                                }
                            }
                        }
                    });
        }

        /**
         * Groups reads by keys - keys are tuples of (library, contig, position, orientation).
         */
        PTransform<PCollection<? extends Read>, PCollection<KV<String, Read>>> makeKeysForFragments(final PCollectionView<SAMFileHeader> headerPcolView) {
            return ParDo
                    .named("make keys for reads")
                    .withSideInputs(headerPcolView)
                    .of(new DoFn<Read, KV<String, Read>>() {
                        private static final long serialVersionUID = 1L;
                        @Override
                        public void processElement(final ProcessContext context) throws Exception {
                            final Read record = context.element();
                            record.setDuplicateFragment(false);
                            final SAMFileHeader h = context.sideInput(headerPcolView);
                            final String key = MarkDuplicatesReadsKey.keyForFragment(h, record);
                            final KV<String, Read> kv = KV.of(key, record);
                            context.output(kv);
                        }
                    });
        }
    }

    /**
     * How to assign a score to the read in MarkDuplicates (so that we pick the best one to be the non-duplicate).
     */
    //Note: copied from htsjdk.samtools.DuplicateScoringStrategy
    private static int score(final Read record) {
        if (record == null) {
            return 0;
        } else {
            return record.getAlignedQuality().stream().mapToInt(Integer::intValue).filter(i -> i >= MIN_BASE_QUAL).sum();
        }
    }

    /**
     * Returns the list of intervals from the given sequence dictionary.
     */
    private static List<SimpleInterval> getAllIntervalsForReference(final SAMSequenceDictionary sequenceDictionary) {
    return GenomeLocSortedSet.createSetFromSequenceDictionary(sequenceDictionary)
            .stream()
            .map(SimpleInterval::new)
            .collect(Collectors.toList());
    }

    /**
     * Struct-like class to store information about the paired reads for mark duplicates.
     */
    private static final class PairedEnds extends GenericJson{
        private Read first, second;

        PairedEnds(final Read first) {
            this.first = first;
        }

        public static PairedEnds of(final Read first) {
            return new PairedEnds(first);
        }

        public PairedEnds and(final Read second) {
            if (second != null && GenomicsReadUtils.unclippedCoordinate(this.first) > GenomicsReadUtils.unclippedCoordinate(second)) {
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

        public Read first() {
            return first;
        }

        public Read second() {
            return second;
        }

        public int score() {
            return MarkDuplicatesDataflow.score(first) + MarkDuplicatesDataflow.score(second);
        }

        @Override
        public PairedEnds clone() {
            return (PairedEnds) super.clone();
        }
    }

    /**
     * Special coder for the PairedEnds class.
     */
    private static final class PairedEndsCoder extends AtomicCoder<PairedEnds> {
        private static final long serialVersionUID = 1L;
        private final KvCoder<Read, Read> delegate =
                KvCoder.of(GenericJsonCoder.of(Read.class), GenericJsonCoder.of(Read.class));

        public static PairedEndsCoder of() {
            return new PairedEndsCoder();
        }

        @Override
        public void encode(final PairedEnds pair, final OutputStream stream, final Context context) throws IOException {
            final Read first = pair.first();
            if (first == null) {
                throw new IOException("nothing to encode");
            }
            Read second = null;
            if (pair.second() != null) {
                second = pair.second();
            } else {
                second = new Read();
                boolean empty = second.isEmpty();//AK: cannot be null because decoding fails
            }
            delegate.encode(KV.of(first, second), stream, context);
        }

        @Override
        public PairedEnds decode(final InputStream stream, final Context context) throws IOException {
            final KV<Read, Read> pair = delegate.decode(stream, context);
            final PairedEnds ends = PairedEnds.of(pair.getKey());
            if (!pair.getValue().isEmpty()) {
                ends.and(pair.getValue());
            }
            return ends;
        }

        @Override
        public void verifyDeterministic() throws NonDeterministicException {
            delegate.verifyDeterministic();
        }
    }

    private static final class CoordinateOrder implements Comparator<Read>, Serializable {
        private static final long serialVersionUID = 1l;
        private final SAMFileHeader header;

        public CoordinateOrder(final SAMFileHeader header) {
            this.header = header;
        }

        @Override
        public int compare(final Read lhs, final Read rhs) {
            if (rhs == lhs) return 0; //shortcut

            final int res1 = Integer.compare(MarkDuplicatesReadsKey.index(header, getReferenceName(lhs)), MarkDuplicatesReadsKey.index(header, getReferenceName(rhs)));
            if (res1 != 0) return res1;

            final int res2 = Long.compare(getAlignmentStart(lhs), getAlignmentStart(rhs));
            if (res2 != 0) return res2;

            final int res3 = Boolean.compare(lhs.getDuplicateFragment(), rhs.getDuplicateFragment());
            if (res3 != 0) return res3;

            final int res4 = Boolean.compare(lhs.getFailedVendorQualityChecks(), rhs.getFailedVendorQualityChecks());
            if (res4 != 0) return res4;

            final int res5 = Integer.compare(lhs.getNumberReads(), rhs.getNumberReads());
            if (res5 != 0) return res5;

            final int res6 = Boolean.compare(lhs.getProperPlacement(), rhs.getProperPlacement());
            if (res6 != 0) return res6;

            final int res7 = Integer.compare(lhs.getReadNumber(), rhs.getReadNumber());
            if (res7 != 0) return res7;

            final int res8 = Boolean.compare(lhs.getSecondaryAlignment(), rhs.getSecondaryAlignment());
            if (res8 != 0) return res8;

            final int res9 = Boolean.compare(lhs.getSupplementaryAlignment(), rhs.getSupplementaryAlignment());
            if (res9 != 0) return res9;

            final int res10 = Integer.compare(GenomicsReadUtils.getMappingQuality(lhs), GenomicsReadUtils.getMappingQuality(rhs));
            if (res10 != 0) return res10;

            final int res11 = Integer.compare(MarkDuplicatesReadsKey.index(header, GenomicsReadUtils.getMateReferenceName(lhs)), MarkDuplicatesReadsKey.index(header, GenomicsReadUtils.getMateReferenceName(rhs)));
            if (res11 != 0) return res11;

            final int res12 = Long.compare(GenomicsReadUtils.getMateAlignmentStart(lhs), GenomicsReadUtils.getMateAlignmentStart(rhs));
            return res12 ;
        }
    }
}