package org.broadinstitute.hellbender.tools.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.*;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.transforms.Partition.PartitionFn;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionList;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.common.collect.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.utils.dataflow.SmallBamWriter;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.tools.dataflow.transforms.MarkDuplicatesReadsKey;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.DuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.LibraryIdGenerator;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

@CommandLineProgramProperties(
        summary ="Marks duplicates on dataflow",
        oneLineSummary ="Mark Duplicates",
        programGroup = DataFlowProgramGroup.class)
public final class MarkDuplicatesDataflow extends DataflowCommandLineProgram {
    private static final long serialVersionUID = 1L;

    //Bases below this quality will not be included in picking the best read from a set of duplicates.
    private static final int MIN_BASE_QUAL = 15;

    // Used to set an attribute on the GATKRead marking this read as an optical duplicate.
    private static final String OPTICAL_DUPLICATE_ATTRIBUTE_NAME = "OD";

    @Argument(doc="output BAM file", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "uri for the input bam, either a local file path or a gs:// bucket path",
              shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
              optional = false)
    protected String bam;

    @Argument(doc = "File to write duplication metrics to.", optional=true,
              shortName = "M", fullName = "METRICS_FILE")
    protected File metricsFile;

    @Argument(doc = "Regular expression that can be used to parse read names in the incoming SAM file. Read names are " +
            "parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used " +
            "to estimate the rate of optical duplication in order to give a more accurate estimated library size. " +
            "Set this option to null to disable optical duplicate detection. " +
            "The regular expression should contain three capture groups for the three variables, in order. " +
            "It must match the entire read name. " +
            "Note that if the default regex is specified, a regex match is not actually done, but instead the read name " +
            " is split on colon character. " +
            "For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. " +
            "For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values.",
            optional = true)
    public String READ_NAME_REGEX = OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX;

    @Argument(doc = "The maximum offset between two duplicate clusters in order to consider them optical duplicates. This " +
            "should usually be set to some fairly small number (e.g. 5-10 pixels) unless using later versions of the " +
            "Illumina pipeline that multiply pixel values by 10, in which case 50-100 is more normal.")
    public int OPTICAL_DUPLICATE_PIXEL_DISTANCE = OpticalDuplicateFinder.DEFAULT_OPTICAL_DUPLICATE_DISTANCE;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Override
    protected void setupPipeline(final Pipeline pipeline) {
        final ReadsDataflowSource readsSource = new ReadsDataflowSource(bam, pipeline);
        final SAMFileHeader header = readsSource.getHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(sequenceDictionary):
                getAllIntervalsForReference(sequenceDictionary);

        final PCollectionView<SAMFileHeader> headerPcolView = pipeline.apply(Create.of(header)).apply(View.<SAMFileHeader>asSingleton());

        final PCollection<GATKRead> preads = readsSource.getReadPCollection(intervals);

        final OpticalDuplicateFinder finder = READ_NAME_REGEX != null ?
            new OpticalDuplicateFinder(OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX, OPTICAL_DUPLICATE_PIXEL_DISTANCE, null) : null;

        final PCollection<GATKRead> results = preads.apply(new MarkDuplicatesDataflowTransform(headerPcolView, finder));

        // TODO: support writing large output files (need a sharded BAM writer)
        SmallBamWriter.writeToFile(pipeline, results, header, outputFile);

        if (metricsFile != null) {
            final PCollection<KV<String,DuplicationMetrics>> metrics = results.apply(new GenerateMetricsTransform(headerPcolView));
            writeMetricsToFile(pipeline, metrics, header, metricsFile);
        }
    }

    /**
     * The main trasform for MarkDuplicates. Takes reads, returned reads with some reads marked as duplicates.
     */
    private static final class MarkDuplicatesDataflowTransform extends PTransform<PCollection<GATKRead>, PCollection<GATKRead>> {

        private static final long serialVersionUID = 1l;

        //private vs non-primary alignments
        private enum ReadsPartition {
            PRIMARY, NOT_PRIMARY
        }

        /**
         *  The header for all reads processed.
        */
        private final PCollectionView<SAMFileHeader> header;

        /**
         * Optical duplicate finder finds the physical location of the reads and determines if they are close enough to
         * be considered optical duplicates.
         */
        private final OpticalDuplicateFinder opticalDuplicateFinder;

        public MarkDuplicatesDataflowTransform( final PCollectionView<SAMFileHeader> header,
            final OpticalDuplicateFinder opticalDuplicateFinder) {
            this.header = header;
            this.opticalDuplicateFinder = opticalDuplicateFinder;
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
            return preads.apply(Partition.of(2, new PartitionFn<GATKRead>() {
                private static final long serialVersionUID = 1l;

                @Override
                public int partitionFor(final GATKRead read, final int n) {
                    return read.isSecondaryAlignment() || read.isSupplementaryAlignment() || read.isUnmapped() ? ReadsPartition.NOT_PRIMARY.ordinal() : ReadsPartition.PRIMARY.ordinal();
                }
            }));
        }

        /**
         * (1) labelReadsByName: label each read with its read group and read name.
         * (2) GroupByKey: group together reads with the same group and name.
         * (3) labelPairedEndsWithAlignmentInfo:
         *   (a) Sort each group of reads (see LexicographicSort below).
         *   (b) Pair consecutive reads into PairedEnds. In most cases there will only be two reads
         *       with the same name. TODO: explain why there might be more.
         *   (c) Label each read with alignment information: Library, reference index,
         *       stranded unclipped start and reverse strand.
         *   (d) Leftover reads are emitted, unmodified, as an unpaired end.
         * (4) GroupByKey: Group PairedEnds that share alignment information. These pairs
         *     are duplicates of each other.
         * (5) markDuplicatePairs:
         *   (a) For each group created by (4), sort the pairs by score and mark all but the
         *       highest scoring as duplicates.
         *   (b) Determine which duplicates are optical duplicates and increase the overall count.
         */
        private PCollection<GATKRead> transformReads(final PCollectionView<SAMFileHeader> headerPcolView, final PCollection<GATKRead> pairs) {
            final PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> labelReadsByName = labelReadsByName(headerPcolView);
            final PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<KV<String, PairedEnds>>> labelPairedEndsWithAlignmentInfo =
                labelPairedEndsWithAlignmentInfo(headerPcolView);
            final PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<GATKRead>> markDuplicatePairedEnds = markDuplicatePairedEnds();

            return pairs
                    .apply(labelReadsByName)
                    .apply(GroupByKey.<String, GATKRead>create())
                    .apply(labelPairedEndsWithAlignmentInfo)
                    .setCoder(KvCoder.of(StringUtf8Coder.of(), new PairedEndsCoder()))
                    .apply(GroupByKey.<String, PairedEnds>create())
                    .apply(markDuplicatePairedEnds);
        }

        private PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> labelReadsByName(final PCollectionView<SAMFileHeader> headerPcolView) {
            return ParDo
                    .named("label reads by name")
                    .withSideInputs(headerPcolView)
                    .of(new DoFn<GATKRead, KV<String, GATKRead>>() {
                        private static final long serialVersionUID = 1l;

                        @Override
                        public void processElement(final ProcessContext context) throws Exception {
                            final GATKRead record = context.element();
                            if (ReadUtils.readHasMappedMate(record)) {
                                final SAMFileHeader h = context.sideInput(headerPcolView);
                                final String key = MarkDuplicatesReadsKey.keyForRead(h, record);
                                final KV<String, GATKRead> kv = KV.of(key, record);
                                context.output(kv);
                            }
                        }
                    });
        }

        private PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<KV<String, PairedEnds>>> labelPairedEndsWithAlignmentInfo(final PCollectionView<SAMFileHeader> headerPcolView) {
            return ParDo
                    .named("label paired ends with alignment info")
                    .withSideInputs(headerPcolView)
                    .of(new DoFn<KV<String, Iterable<GATKRead>>, KV<String, PairedEnds>>() {
                        private static final long serialVersionUID = 1L;
                        @Override
                        public void processElement(final ProcessContext context) throws Exception {
                            final SAMFileHeader header = context.sideInput(headerPcolView);
                            final List<GATKRead> sorted = Lists.newArrayList(context.element().getValue());
                            sorted.sort(new LexicographicOrder(header));
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


        private PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<GATKRead>> markDuplicatePairedEnds() {
            return ParDo
                    .named("mark duplicate paired ends")
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
                            List<PairedEnds> scored = Ordering.natural().reverse().onResultOf((PairedEnds pair) -> pair.score()).sortedCopy(paired.get(true));
                            final PairedEnds best = Iterables.getFirst(scored, null);
                            if (best == null) {
                                return;
                            }
                            final String bestName = best.first().getName();
                            for (PairedEnds pair : scored) {
                                opticalDuplicateFinder.addLocationInformation(pair.first().getName(), pair);
                            }

                            // We do not need to split the list by orientation as the keys for the pairs already
                            // include directionality information and a FR pair would not be grouped with an RF pair.
                            // This function sorts the scored list, so we keep track of the name of the best one,
                            // so that it is not marked as a duplicate.
                            final boolean[] opticalDuplicateFlags = opticalDuplicateFinder.findOpticalDuplicates(scored);

                            //Mark everyone who's not best as a duplicate
                            for (int i = 0; i < scored.size(); ++i) {
                                PairedEnds pair = scored.get(i);
                                String pairName = pair.first().getName();

                                GATKRead record = pair.first();
                                if (!pairName.equals(bestName)) {
                                    record.setIsDuplicate(true);
                                }

                                // Note: The read marked as an optical duplicate may not be marked as a duplicate. The reason
                                // for this is that the optical duplicate finder sorts based on a different criterion then
                                // the one above. We really only care about the total number of optical duplicates, so this
                                // is OK.
                                if (opticalDuplicateFlags[i]) {
                                    record.setAttribute(OPTICAL_DUPLICATE_ATTRIBUTE_NAME, 1);
                                }
                                context.output(record);

                                record = pair.second();
                                if (!pairName.equals(bestName)) {
                                    record.setIsDuplicate(true);
                                }
                                if (opticalDuplicateFlags[i]) {
                                    record.setAttribute(OPTICAL_DUPLICATE_ATTRIBUTE_NAME, 1);
                                }
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
         * Groups reads by keys - keys are tuples of (library, contig, position, orientation).
         */
        PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> makeKeysForFragments(final PCollectionView<SAMFileHeader> headerPcolView) {
            return ParDo
                    .named("make keys for fragments")
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

    // This transform takes marked reads and generates the metrics from those tests.
    private static final class GenerateMetricsTransform extends PTransform<PCollection<GATKRead>, PCollection<KV<String, DuplicationMetrics>>> {
        private static final long serialVersionUID = 1l;

        /**
         *  The header for all reads processed.
        */
        private final PCollectionView<SAMFileHeader> header;

        public GenerateMetricsTransform(final PCollectionView<SAMFileHeader> header) {
            this.header = header;
        }

        @Override
        public PCollection<KV<String, DuplicationMetrics>> apply(final PCollection<GATKRead> preads) {
            return preads
                .apply(generateMetricsByLibrary())
                .setCoder(KvCoder.of(StringUtf8Coder.of(), new DuplicationMetricsCoder()))
                // Combine metrics for reads from the same library.
                .apply(Combine.<String, DuplicationMetrics>perKey(new CombineMetricsFn()))
                // For each library, finalize the metrics by computing derived metrics and dividing paired counters
                // by 2.
                .apply(finalizeMetrics());
        }

        private PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, DuplicationMetrics>>> generateMetricsByLibrary() {
            return ParDo
                .named("generate metrics by library")
                .withSideInputs(header)
                .of(new DoFn<GATKRead, KV<String, DuplicationMetrics>>() {
                        private static final long serialVersionUID = 1l;

                        @Override
                        public void processElement(final ProcessContext context) throws Exception {
                            final GATKRead record = context.element();
                            final SAMFileHeader h = context.sideInput(header);
                            final String library = LibraryIdGenerator.getLibraryName(h, record.getReadGroup());
                            DuplicationMetrics metrics = new DuplicationMetrics();
                            metrics.LIBRARY = library;
                            // These need to be initialized so that the DuplicationMetrics can be serialized safely.
                            metrics.PERCENT_DUPLICATION = 0.0;
                            metrics.ESTIMATED_LIBRARY_SIZE = 0L;
                            // We have already filtered out secondary and supplementary reads in the partition transform.
                            if (record.isUnmapped()) {
                                ++metrics.UNMAPPED_READS;
                            } else if (!record.isPaired() || record.mateIsUnmapped()) {
                                ++metrics.UNPAIRED_READS_EXAMINED;
                            } else {
                                ++metrics.READ_PAIRS_EXAMINED;
                            }

                            if (record.isDuplicate()) {
                                if (!record.isPaired() || record.mateIsUnmapped()) {
                                    ++metrics.UNPAIRED_READ_DUPLICATES;
                                } else {
                                    ++metrics.READ_PAIR_DUPLICATES;
                                }
                            }
                            if (record.hasAttribute(OPTICAL_DUPLICATE_ATTRIBUTE_NAME) &&
                                record.getAttributeAsInteger(OPTICAL_DUPLICATE_ATTRIBUTE_NAME) == 1) {
                                ++metrics.READ_PAIR_OPTICAL_DUPLICATES;
                            }
                            final KV<String, DuplicationMetrics> kv = KV.of(library, metrics);
                            context.output(kv);
                        }
                    });
        }


        public static class CombineMetricsFn implements SerializableFunction<Iterable<DuplicationMetrics>, DuplicationMetrics> {
            private static final long serialVersionUID = 1l;

            @Override
            public DuplicationMetrics apply(Iterable<DuplicationMetrics> input) {
                DuplicationMetrics metricsSum = new DuplicationMetrics();
                int num = 0;
                for (final DuplicationMetrics metrics : input) {
                    if (metricsSum.LIBRARY == null) {
                        metricsSum.LIBRARY = metrics.LIBRARY;
                    }
                    // This should never happen, as we grouped by key using library as the key.
                    if (!metricsSum.LIBRARY.equals(metrics.LIBRARY)) {
                        return null;
                    }
                    metricsSum.UNMAPPED_READS += metrics.UNMAPPED_READS;
                    metricsSum.UNPAIRED_READS_EXAMINED += metrics.UNPAIRED_READS_EXAMINED;
                    metricsSum.READ_PAIRS_EXAMINED += metrics.READ_PAIRS_EXAMINED;
                    metricsSum.UNPAIRED_READ_DUPLICATES += metrics.UNPAIRED_READ_DUPLICATES;
                    metricsSum.READ_PAIR_DUPLICATES += metrics.READ_PAIR_DUPLICATES;
                    metricsSum.READ_PAIR_OPTICAL_DUPLICATES += metrics.READ_PAIR_OPTICAL_DUPLICATES;
                }
                metricsSum.PERCENT_DUPLICATION = 0.0;
                metricsSum.ESTIMATED_LIBRARY_SIZE = 0L;
                return metricsSum;
            }
        }

        private PTransform<PCollection<? extends KV<String, DuplicationMetrics>>, PCollection<KV<String, DuplicationMetrics>>> finalizeMetrics() {
            return ParDo
                .named("finalize metrics")
                .of(new DoFn<KV<String, DuplicationMetrics>, KV<String, DuplicationMetrics>>() {
                        private static final long serialVersionUID = 1l;

                        @Override
                        public void processElement(final ProcessContext context) throws Exception {
                            DuplicationMetrics metrics = context.element().getValue();
                            metrics.READ_PAIRS_EXAMINED = metrics.READ_PAIRS_EXAMINED / 2;
                            metrics.READ_PAIR_DUPLICATES = metrics.READ_PAIR_DUPLICATES / 2;
                            metrics.READ_PAIR_OPTICAL_DUPLICATES = metrics.READ_PAIR_OPTICAL_DUPLICATES / 2;
                            metrics.calculateDerivedMetrics();
                            if (metrics.ESTIMATED_LIBRARY_SIZE == null) {
                                metrics.ESTIMATED_LIBRARY_SIZE = 0L;
                            }
                            final KV<String, DuplicationMetrics> kv = KV.of(context.element().getKey(), metrics);
                            context.output(kv);
                        }
                    });
        }
    }


    /**
     * Special coder for the DuplicationMetrics class.
     */
    private static final class DuplicationMetricsCoder extends CustomCoder<DuplicationMetrics> {
        private static final long serialVersionUID = 1L;

        @Override
        public void encode( DuplicationMetrics value, OutputStream outStream, Context context ) throws CoderException, IOException {
            if ( value == null ) {
                throw new IOException("nothing to encode");
            }
            SerializableCoder.of(String.class).encode(value.LIBRARY, outStream, context);
            SerializableCoder.of(Long.class).encode(value.UNMAPPED_READS, outStream, context);
            SerializableCoder.of(Long.class).encode(value.UNPAIRED_READS_EXAMINED, outStream, context);
            SerializableCoder.of(Long.class).encode(value.READ_PAIRS_EXAMINED, outStream, context);
            SerializableCoder.of(Long.class).encode(value.UNPAIRED_READ_DUPLICATES, outStream, context);
            SerializableCoder.of(Long.class).encode(value.READ_PAIR_DUPLICATES, outStream, context);
            SerializableCoder.of(Long.class).encode(value.READ_PAIR_OPTICAL_DUPLICATES, outStream, context);
            SerializableCoder.of(Double.class).encode(value.PERCENT_DUPLICATION, outStream, context);
            SerializableCoder.of(Long.class).encode(value.ESTIMATED_LIBRARY_SIZE, outStream, context);
        }

        @Override
        public DuplicationMetrics decode( InputStream inStream, Context context ) throws CoderException, IOException {
            DuplicationMetrics metrics = new DuplicationMetrics();
            metrics.LIBRARY = SerializableCoder.of(String.class).decode(inStream, context);
            metrics.UNMAPPED_READS = SerializableCoder.of(Long.class).decode(inStream, context);
            metrics.UNPAIRED_READS_EXAMINED = SerializableCoder.of(Long.class).decode(inStream, context);
            metrics.READ_PAIRS_EXAMINED = SerializableCoder.of(Long.class).decode(inStream, context);
            metrics.UNPAIRED_READ_DUPLICATES = SerializableCoder.of(Long.class).decode(inStream, context);
            metrics.READ_PAIR_DUPLICATES = SerializableCoder.of(Long.class).decode(inStream, context);
            metrics.READ_PAIR_OPTICAL_DUPLICATES = SerializableCoder.of(Long.class).decode(inStream, context);
            metrics.PERCENT_DUPLICATION = SerializableCoder.of(Double.class).decode(inStream, context);
            metrics.ESTIMATED_LIBRARY_SIZE = SerializableCoder.of(Long.class).decode(inStream, context);
            return metrics;
        }
    }

    private static void writeMetricsToFile(Pipeline pipeline, final PCollection<KV<String,DuplicationMetrics>> metrics,
                                           final SAMFileHeader header, final File dest) {
        PCollectionView<Map<String,Iterable<DuplicationMetrics>>> mapView = metrics
            .apply(View.asMap());

        PCollection<File> dummy = pipeline.apply(Create.<File>of(dest).setName("output metrics file"));

        dummy.apply(ParDo.named("save metrics file")
                    .withSideInputs(mapView)
                    .of(new WriteToMetricsFile(header, mapView))
        );

    }

    private static class WriteToMetricsFile extends DoFn<File,Void> implements Serializable {
        private static final long serialVersionUID = 1L;
        private final SAMFileHeader header;
        private final PCollectionView<Map<String,Iterable<DuplicationMetrics>>> mapView;

        public WriteToMetricsFile(SAMFileHeader header, PCollectionView<Map<String,Iterable<DuplicationMetrics>>> mapView) {
            this.header = header;
            this.mapView = mapView;
        }

        @Override
        public void processElement(ProcessContext c) throws Exception {
            File dest = c.element();
            Map<String,Iterable<DuplicationMetrics>> metrics = c.sideInput(mapView);
            final MetricsFile<DuplicationMetrics,Double> file = new MetricsFile<>();
            for (final Map.Entry<String,Iterable<DuplicationMetrics>> entry : metrics.entrySet()) {
                for (final DuplicationMetrics m : entry.getValue()) {
                    file.addMetric(m);
                }
            }
            // Should we add a histogram?
            file.write(dest);
        }
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
    private static final class PairedEnds implements OpticalDuplicateFinder.PhysicalLocation {
        private GATKRead first, second;

        // Information used to detect optical dupes
        public short readGroup = -1;
        public short tile = -1;
        public short x = -1, y = -1;
        public short libraryId = -1;

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
            return MarkDuplicatesDataflow.score(first) + MarkDuplicatesDataflow.score(second);
        }

        @Override
        public short getReadGroup() { return this.readGroup; }

        @Override
        public void setReadGroup(final short readGroup) { this.readGroup = readGroup; }

        @Override
        public short getTile() { return this.tile; }

        @Override
        public void setTile(final short tile) { this.tile = tile; }

        @Override
        public short getX() { return this.x; }

        @Override
        public void setX(final short x) { this.x = x; }

        @Override
        public short getY() { return this.y; }

        @Override
        public void setY(final short y) { this.y = y; }

        @Override
        public short getLibraryId() { return this.libraryId; }

        @Override
        public void setLibraryId(final short libraryId) { this.libraryId = libraryId; }
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
    // TODO: explain how this ordering works.
    private static final class LexicographicOrder implements Comparator<GATKRead>, Serializable {
        private static final long serialVersionUID = 1l;
        private final SAMFileHeader header;

        public LexicographicOrder(final SAMFileHeader header) {
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
}
