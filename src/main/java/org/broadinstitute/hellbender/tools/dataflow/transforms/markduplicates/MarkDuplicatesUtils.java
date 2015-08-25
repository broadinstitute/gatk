package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.coders.StringUtf8Coder;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.Filter;
import com.google.cloud.dataflow.sdk.transforms.GroupByKey;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.common.collect.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.dataflow.coders.PairedEndsCoder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.DuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.LibraryIdGenerator;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

import java.io.File;
import java.io.InputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Utility classes and functions for Mark Duplicates.
 */
final class MarkDuplicatesUtils {
    //Bases below this quality will not be included in picking the best read from a set of duplicates.
    public static final int MIN_BASE_QUAL = 15;

    // Used to set an attribute on the GATKRead marking this read as an optical duplicate.
    private static final String OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME = "OD";

  /**
   * (1) keyReadsByName: label each read with its read group and read name.
   * (2) GroupByKey: group together reads with the same group and name.
   * (3) keyPairedEndsWithAlignmentInfo:
   *   (a) Sort each group of reads (see GATKOrder below).
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
  static PCollection<GATKRead> transformReads(final PCollectionView<SAMFileHeader> headerPcolView, final PCollectionView<OpticalDuplicateFinder> finderPcolView, final PCollection<GATKRead> pairs) {

      final SerializableFunction<GATKRead, Boolean> hasMappedMate = new SerializableFunction<GATKRead, Boolean>() {
          @Override
          public Boolean apply(GATKRead read) {
              return ReadUtils.readHasMappedMate(read);
          }
      } ;
      final PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> keyReadsByName = keyReadsByName(headerPcolView);
      final PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<KV<String, PairedEnds>>> keyPairedEndsWithAlignmentInfo =
          keyPairedEndsWithAlignmentInfo(headerPcolView);
      final PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<GATKRead>> markDuplicatePairedEnds = markDuplicatePairedEnds(finderPcolView);

      return pairs
        .apply(Filter.by(hasMappedMate))
        .apply(keyReadsByName)
        .apply(GroupByKey.<String, GATKRead>create())
        .apply(keyPairedEndsWithAlignmentInfo)
        .setCoder(KvCoder.of(StringUtf8Coder.of(), new PairedEndsCoder()))
        .apply(GroupByKey.<String, PairedEnds>create())
        .apply(markDuplicatePairedEnds);
  }

    /**
     * Makes keys for read pairs. To be grouped by in the next step.
     */
    static PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> keyReadsByName(final PCollectionView<SAMFileHeader> headerPcolView) {
        BiFunction<SAMFileHeader, GATKRead ,String> func = (BiFunction<SAMFileHeader, GATKRead ,String> & Serializable) MarkDuplicatesUtils::keyForRead;
        return applyKeyFunction(headerPcolView, func, "key reads by name");
    }

    // the specific key
    static String keyForRead(SAMFileHeader header, GATKRead read) {
        return ReadsKey.keyForRead(header, read);
    }

    // given reads that have the same ReadKey (i.e. are at the same location), we output PairedEnds.
    static ArrayList<PairedEnds> pairTheReads(SAMFileHeader header, Iterable<GATKRead> keyedReads) {
        ArrayList<PairedEnds> ret = new ArrayList<>();
        final List<GATKRead> sorted = Lists.newArrayList(keyedReads);
        sorted.sort(new GATKOrder(header));
        PairedEnds pair = null;
        //Records are sorted, we iterate over them and pair them up.
        for (final GATKRead record : sorted) {
            if (pair == null) {                                //first in pair
                pair = PairedEnds.of(record);
            } else {                                           //second in pair
                pair.and(record);
                ret.add(pair);
                pair = null;                                   //back to first
            }
        }
        if (pair != null) {                                    //left over read
            ret.add(pair);
        }
        return ret;
    }

    // PairedEnd -> its key
    static String keyForPair(SAMFileHeader header, PairedEnds pair) {
        return pair.key(header);
    }

    // Dataflowese for "pairTheReads and then keyForPair"
    static PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<KV<String, PairedEnds>>> keyPairedEndsWithAlignmentInfo(final PCollectionView<SAMFileHeader> headerPcolView) {
        return ParDo
                .named("key paired ends with alignment info")
                .withSideInputs(headerPcolView)
                .of(new DoFn<KV<String, Iterable<GATKRead>>, KV<String, PairedEnds>>() {
                    private static final long serialVersionUID = 1L;

                    @Override
                    public void processElement(final ProcessContext context) throws Exception {
                        final SAMFileHeader header = context.sideInput(headerPcolView);
                        for (PairedEnds pair : pairTheReads(header, context.element().getValue())) {
                            context.output(KV.of(keyForPair(header, pair), pair));
                        }
                    }
                });
    }

    static ArrayList<GATKRead> markDuplicatePairedEndsFn(OpticalDuplicateFinder finder, Iterable<PairedEnds> pairs) {
        ArrayList<GATKRead> ret = new ArrayList<>();
        final Iterable<PairedEnds> pairedEndsCopy = Iterables.transform(pairs, PairedEnds::copy);
        final ImmutableListMultimap<Boolean, PairedEnds> paired = Multimaps.index(pairedEndsCopy, pair -> pair.second() != null);

        // As in Picard, unpaired ends left alone.
        for (final PairedEnds pair : paired.get(false)) {
            ret.add(pair.first());
        }

        // order by score
        List<PairedEnds> scored = Ordering.natural().reverse().onResultOf((PairedEnds pair) -> pair.score()).sortedCopy(paired.get(true));

        final PairedEnds best = Iterables.getFirst(scored, null);
        if (best == null) {
            return ret;
        }

        // Mark everyone who's not best as a duplicate
        for (final PairedEnds pair : Iterables.skip(scored, 1)) {
            pair.first().setIsDuplicate(true);
            pair.second().setIsDuplicate(true);
        }
        // Now, add location information to the paired ends
        for (final PairedEnds pair : scored) {
            // Both elements in the pair have the same name
            finder.addLocationInformation(pair.first().getName(), pair);
        }

        // This must happen last, as findOpticalDuplicates mutates the list.
        // We do not need to split the list by orientation as the keys for the pairs already
        // include directionality information and a FR pair would not be grouped with an RF pair.
        final boolean[] opticalDuplicateFlags = finder.findOpticalDuplicates(scored);
        int numOpticalDuplicates = 0;
        for (final boolean b : opticalDuplicateFlags) {
            if (b) {
                numOpticalDuplicates++;
            }
        }
        best.first().setAttribute(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME, numOpticalDuplicates);

        for (final PairedEnds pair : scored) {
            ret.add(pair.first());
            ret.add(pair.second());
        }
        return ret;
    }

    static PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<GATKRead>> markDuplicatePairedEnds(PCollectionView<OpticalDuplicateFinder> finderPcolView) {
        return ParDo
                .named("mark duplicate paired ends")
                .withSideInputs(finderPcolView)
                .of(new DoFn<KV<String, Iterable<PairedEnds>>, GATKRead>() {
                    private static final long serialVersionUID = 1l;

                    @Override
                    public void processElement(final ProcessContext context) throws Exception {
                        final OpticalDuplicateFinder finder = context.sideInput(finderPcolView);
                        for (GATKRead read : markDuplicatePairedEndsFn(finder, context.element().getValue())) {
                            context.output(read);
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
    static PCollection<GATKRead> transformFragments(final PCollectionView<SAMFileHeader> headerPcolView, final PCollection<GATKRead> fragments) {
        final PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> makeKeysForFragments =  makeKeysForFragments(headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<GATKRead>> markGroupedDuplicateFragments = markGroupedDuplicateFragments();
        return fragments
                .apply(markNonDuplicateTransform())
                .apply(makeKeysForFragments)
                .apply(GroupByKey.<String, GATKRead>create())
                .apply(markGroupedDuplicateFragments);//no need to set up coder for Read (uses GenericJsonCoder)
    }

    static ArrayList<GATKRead> transformFragmentsFn(Iterable<GATKRead> input) {
        ArrayList<GATKRead> ret = new ArrayList<>();
        //split reads by paired vs unpaired
        Iterable<GATKRead> readsCopy = Iterables.transform(input, GATKRead::copy);
        final Map<Boolean, List<GATKRead>> byPairing = StreamSupport.stream(readsCopy.spliterator(), false).collect(Collectors.partitioningBy(
            read -> ReadUtils.readHasMappedMate(read)
        ));

        // Note the we emit only fragments from this mapper.
        if (byPairing.get(true).isEmpty()) {
            // There are no paired reads, mark all but the highest scoring fragment as duplicate.
            final List<GATKRead> frags = Ordering.natural().reverse().onResultOf((GATKRead read) -> MarkDuplicatesUtils.score(read)).immutableSortedCopy(byPairing.get(false));
            if (!frags.isEmpty()) {
                ret.add(frags.get(0));                        //highest score - just emit
                for (final GATKRead record : Iterables.skip(frags, 1)) {  //lower   scores - mark as dups and emit
                    record.setIsDuplicate(true);
                    ret.add(record);
                }
            }
        } else {
            // There are paired ends so we mark all fragments as duplicates.
            for (final GATKRead record : byPairing.get(false)) {
                record.setIsDuplicate(true);
                ret.add(record);
            }
        }
        return ret;
    }

    static PTransform<PCollection<? extends KV<String, Iterable<GATKRead>>>, PCollection<GATKRead>> markGroupedDuplicateFragments() {
        return ParDo.named("mark dups")
                .of(new DoFn<KV<String, Iterable<GATKRead>>, GATKRead>() {
                    private static final long serialVersionUID = 1l;

                    @Override
                    public void processElement(final ProcessContext context) throws Exception {
                        for (GATKRead r : transformFragmentsFn(context.element().getValue())) {
                            context.output(r);
                        }
                    }
                });
    }



    // the specific key
    static String keyForFragment(SAMFileHeader header, GATKRead read) {
        return ReadsKey.keyForFragment(header, read);
    }


    /**
     * Groups reads by keys - keys are tuples of (library, contig, position, orientation).
     */
    static PTransform<PCollection<? extends GATKRead>, PCollection<KV<String, GATKRead>>> makeKeysForFragments(final PCollectionView<SAMFileHeader> headerPcolView) {
        BiFunction<SAMFileHeader, GATKRead ,String> func = (BiFunction<SAMFileHeader, GATKRead ,String> & Serializable) MarkDuplicatesUtils::keyForFragment;
        return applyKeyFunction(headerPcolView, func, "key fragments");
    }

    public static GATKRead markNonDuplicate(GATKRead read) {
        final GATKRead record = read.copy();
        record.setIsDuplicate(false);
        return record;
    }

    static PTransform<PCollection<? extends GATKRead>, PCollection<GATKRead>> markNonDuplicateTransform() {
        return ParDo
            .named("mark non-duplicate")
            .of(new DoFn<GATKRead, GATKRead>() {
                private static final long serialVersionUID = 1L;

                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    context.output(markNonDuplicate(context.element()));
                }
            });
    }

    // This transform takes marked reads and generates the metrics from those tests.
    public static class GenerateMetricsTransform extends PTransform<PCollection<GATKRead>, PCollection<KV<String, DuplicationMetrics>>> {
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
                .apply(Filter.by(
                    new SerializableFunction<GATKRead, Boolean>() {
                      private static final long serialVersionUID = 0;

                      public Boolean apply(GATKRead element) {
                        return !element.isSecondaryAlignment() && !element.isSupplementaryAlignment();
                      }
                    }))
                .apply(generateMetricsByLibrary())
                .setCoder(KvCoder.of(StringUtf8Coder.of(), SerializableCoder.of(DuplicationMetrics.class)))
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
                            if (record.isSecondaryAlignment() || record.isSupplementaryAlignment()) {
                              return;
                            }
                            final SAMFileHeader h = context.sideInput(header);
                            final String library = LibraryIdGenerator.getLibraryName(h, record.getReadGroup());
                            DuplicationMetrics metrics = new DuplicationMetrics();
                            metrics.LIBRARY = library;
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
                            if (record.hasAttribute(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME)) {
                              metrics.READ_PAIR_OPTICAL_DUPLICATES +=
                                  record.getAttributeAsInteger(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME);
                            }
                            final KV<String, DuplicationMetrics> kv = KV.of(library, metrics);
                            context.output(kv);
                        }
                    });
        }


        public class CombineMetricsFn implements SerializableFunction<Iterable<DuplicationMetrics>, DuplicationMetrics> {
            private static final long serialVersionUID = 1l;

            @Override
            public DuplicationMetrics apply(Iterable<DuplicationMetrics> input) {
                DuplicationMetrics metricsSum = new DuplicationMetrics();
                for (final DuplicationMetrics metrics : input) {
                    if (metricsSum.LIBRARY == null) {
                        metricsSum.LIBRARY = metrics.LIBRARY;
                    }
                    // This should never happen, as we grouped by key using library as the key.
                    if (!metricsSum.LIBRARY.equals(metrics.LIBRARY)) {
                      throw new GATKException("Two different libraries encountered while summing metrics: " + metricsSum.LIBRARY
                          + " and " + metrics.LIBRARY);
                    }
                    metricsSum.UNMAPPED_READS += metrics.UNMAPPED_READS;
                    metricsSum.UNPAIRED_READS_EXAMINED += metrics.UNPAIRED_READS_EXAMINED;
                    metricsSum.READ_PAIRS_EXAMINED += metrics.READ_PAIRS_EXAMINED;
                    metricsSum.UNPAIRED_READ_DUPLICATES += metrics.UNPAIRED_READ_DUPLICATES;
                    metricsSum.READ_PAIR_DUPLICATES += metrics.READ_PAIR_DUPLICATES;
                    metricsSum.READ_PAIR_OPTICAL_DUPLICATES += metrics.READ_PAIR_OPTICAL_DUPLICATES;
                }
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
                            DuplicationMetrics metrics = context.element().getValue().copy();
                            // Divide these by 2 because they are counted for each read
                            // when they should be counted by pair.
                            metrics.READ_PAIRS_EXAMINED = metrics.READ_PAIRS_EXAMINED / 2;
                            metrics.READ_PAIR_DUPLICATES = metrics.READ_PAIR_DUPLICATES / 2;

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


    public static void writeMetricsToFile(Pipeline pipeline, final PCollection<KV<String,DuplicationMetrics>> metrics,
        final SAMFileHeader header, final File dest) {
        PCollectionView<Map<String,Iterable<DuplicationMetrics>>> mapView = metrics
            .apply(View.asMultimap());

        PCollection<File> dummy = pipeline.apply(Create.<File>of(dest));

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
    static int score(final GATKRead record) {
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

    // ------------------------------------------------------------------------------------------------------------
    // Dataflow helper functions

    // key function -> PTransform that puts this key in front.
    static <T> PTransform<PCollection<? extends T>, PCollection<KV<String, T>>> applyKeyFunction(final PCollectionView<SAMFileHeader> headerPcolView, final BiFunction<SAMFileHeader, T ,String> keyFunc, final String name) {
        final BiFunction<SAMFileHeader, T ,String> serializableKeyFunc = (BiFunction<SAMFileHeader, T ,String> & Serializable) keyFunc;
        return ParDo
            .named(name)
            .withSideInputs(headerPcolView)
            .of(new DoFn<T, KV<String, T>>() {
                private static final long serialVersionUID = 1l;

                @Override
                public void processElement(final ProcessContext context) throws Exception {
                    final T record = context.element();
                    final SAMFileHeader h = context.sideInput(headerPcolView);
                    String key = serializableKeyFunc.apply(h, record);
                    if (null != key) context.output(KV.of(key, record));
                }
            });
    }

    // ------------------------------------------------------------------------------------------------------------

    /**
     * GATKRead comparator that compares based on mapping position followed by SAM flags.
     */
    final static class GATKOrder implements Comparator<GATKRead>, Serializable {
        private static final long serialVersionUID = 1l;
        private final SAMFileHeader header;
        // TODO: Unify with other comparators in the codebase

        public GATKOrder(final SAMFileHeader header) {
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
            return res12;
        }
    }
}
