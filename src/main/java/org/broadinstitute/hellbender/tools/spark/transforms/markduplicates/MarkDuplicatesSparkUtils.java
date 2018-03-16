package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.serializers.FieldSerializer;
import com.google.common.collect.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.markduplicates.*;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Utility classes and functions for Mark Duplicates.
 */
public class MarkDuplicatesSparkUtils {

    // Used to set an attribute on the GATKRead marking this read as an optical duplicate.
    public static final String OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME = "OD";
    private static final Comparator<PairedEnds> PAIRED_ENDS_SCORE_COMPARATOR = Comparator.comparing(PairedEnds::getScore)
            .thenComparing(PairedEndsCoordinateComparator.INSTANCE.reversed());


    @DefaultSerializer(FieldSerializer.class)
    public static class IndexPair<T>{
        private final T value;
        private final int index;

        public T getValue() {
            return value;
        }

        public int getIndex() {
            return index;
        }

        private IndexPair(T value, int index) {
            this.value = value;
            this.index = index;
        }
    }

    /**
     * (0) filter: remove unpaired reads and reads with an unmapped mate.
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
    static JavaPairRDD<IndexPair<String>, Integer> transformToDuplicateNames(final SAMFileHeader header, final MarkDuplicatesScoringStrategy scoringStrategy, final OpticalDuplicateFinder finder, final JavaRDD<GATKRead>  reads, final int numReducers) {
        //remove reads that are unmapped and have unmapped or no mates
        // we treat these specially and don't mark them as duplicates
        final JavaRDD<GATKRead> mappedReads = reads.filter(ReadFilterLibrary.MAPPED::test);

        final JavaPairRDD<String, Iterable<IndexPair<GATKRead>>> keyedReads = getReadsGroupedByName(header, mappedReads, numReducers);

        // Place all the reads into a single RDD separated
        final JavaPairRDD<Integer, PairedEnds> pairedEnds = keyedReads.flatMapToPair(keyedRead -> {
            final List<Tuple2<Integer, PairedEnds>> out = Lists.newArrayList();

            final Iterator<IndexPair<GATKRead>> sorted = Utils.stream(keyedRead._2())
                    ////// Making The Fragments //////
                    // Make a PairedEnd object with no second read for each fragment (and an empty one for each paired read)
                    .peek(readWithIndex -> {
                        final GATKRead read = readWithIndex.getValue();
                        PairedEnds fragment = (ReadUtils.readHasMappedMate(read)) ?
                                PairedEnds.placeHolder(read, header, readWithIndex.getIndex()) :
                                PairedEnds.newFragment(read, header, readWithIndex.getIndex(), scoringStrategy);

                        out.add(new Tuple2<>(fragment.keyForFragment(header), fragment));
                    })
                    .filter(readWithIndex -> ReadUtils.readHasMappedMate(readWithIndex.getValue()))
                    .sorted(new GATKOrder(header))
                    .iterator();

            ////// Making The Paired Reads //////
            // Write each paired read with a mapped mate as a pair
            while (sorted.hasNext()) {
                final IndexPair<GATKRead> first = sorted.next();

                final PairedEnds pair;
                if (sorted.hasNext()) {
                    final IndexPair<GATKRead> second = sorted.next();
                    pair = PairedEnds.newPair(first.getValue(), second.getValue(), header, second.getIndex(), scoringStrategy);
                } else {
                    pair = PairedEnds.newPairWithMissingSecond(first.getValue(), header, first.getIndex(), scoringStrategy);
                }
                out.add(new Tuple2<>(pair.key(header), pair));
            }
            return out.iterator();
        });

        final JavaPairRDD<Integer, Iterable<PairedEnds>> keyedPairs = pairedEnds.groupByKey(); //TODO make this a proper aggregate by key

        return markPairedEnds(keyedPairs, finder);
    }

    //todo use this instead of keeping all unmapped reads as non-duplicate
    public static boolean readAndMateAreUnmapped(GATKRead read) {
        return read.isUnmapped() && (!read.isPaired() || read.mateIsUnmapped());
    }

    private static JavaPairRDD<String, Iterable<IndexPair<GATKRead>>> getReadsGroupedByName(SAMFileHeader header, JavaRDD<GATKRead> reads, int numReducers) {

        final JavaPairRDD<String, Iterable<IndexPair<GATKRead>>> keyedReads;
        final JavaRDD<IndexPair<GATKRead>> indexedReads = reads.mapPartitionsWithIndex(
                (index, iter) -> Utils.stream(iter).map(read -> new IndexPair<>(read, index)).iterator(), false);
        if (SAMFileHeader.SortOrder.queryname.equals(header.getSortOrder()) || SAMFileHeader.GroupOrder.query.equals(header.getGroupOrder()) ) {
            // reads are already grouped by name, so perform grouping within the partition (no shuffle)
            keyedReads = spanReadsByKey(indexedReads);
        } else {
            // sort by group and name (incurs a shuffle)
            JavaPairRDD<String, IndexPair<GATKRead>> keyReadPairs = indexedReads.mapToPair(read -> new Tuple2<>(ReadsKey.keyForRead(
                    read.getValue()), read));
            keyedReads = keyReadPairs.groupByKey(numReducers);
        }
        return keyedReads;
    }

    static JavaPairRDD<String, Iterable<IndexPair<GATKRead>>> spanReadsByKey(final JavaRDD<IndexPair<GATKRead>> reads) {
        JavaPairRDD<String, IndexPair<GATKRead>> nameReadPairs = reads.mapToPair(read -> new Tuple2<>(read.getValue().getName(), read));
        return spanByKey(nameReadPairs).flatMapToPair(namedRead -> {
            // for each name, separate reads by key (group name)
            List<Tuple2<String, Iterable<IndexPair<GATKRead>>>> out = Lists.newArrayList();
            ListMultimap<String, IndexPair<GATKRead>> multi = LinkedListMultimap.create();
            for (IndexPair<GATKRead> read : namedRead._2()) {
                multi.put(ReadsKey.keyForRead(read.getValue()), read);
            }
            for (String key : multi.keySet()) {
                // list from Multimap is not serializable by Kryo, so put in a new array list
                out.add(new Tuple2<>(key, Lists.newArrayList(multi.get(key))));
            }
            return out.iterator();
        });
    }

    /**
     * Like <code>groupByKey</code>, but assumes that values are already sorted by key, so no shuffle is needed,
     * which is much faster.
     * @param rdd the input RDD
     * @param <K> type of keys
     * @param <V> type of values
     * @return an RDD where each the values for each key are grouped into an iterable collection
     */
    private static <K, V> JavaPairRDD<K, Iterable<V>> spanByKey(JavaPairRDD<K, V> rdd) {
        return rdd.mapPartitionsToPair(MarkDuplicatesSparkUtils::spanningIterator);
    }

    /**
     * An iterator that groups values having the same key into iterable collections.
     * @param iterator an iterator over key-value pairs
     * @param <K> type of keys
     * @param <V> type of values
     * @return an iterator over pairs of keys and grouped values
     */
    static <K, V> Iterator<Tuple2<K, Iterable<V>>> spanningIterator(Iterator<Tuple2<K, V>> iterator) {
        final PeekingIterator<Tuple2<K, V>> iter = Iterators.peekingIterator(iterator);
        return new AbstractIterator<Tuple2<K, Iterable<V>>>() {
            @Override
            protected Tuple2<K, Iterable<V>> computeNext() {
                K key = null;
                List<V> group = Lists.newArrayList();
                while (iter.hasNext()) {
                    if (key == null) {
                        Tuple2<K, V> next = iter.next();
                        key = next._1();
                        V value = next._2();
                        group.add(value);
                        continue;
                    }
                    K nextKey = iter.peek()._1(); // don't advance...
                    if (nextKey.equals(key)) {
                        group.add(iter.next()._2()); // .. unless the keys match
                    } else {
                        return new Tuple2<>(key, group);
                    }
                }
                if (key != null) {
                    return new Tuple2<>(key, group);
                }
                return endOfData();
            }
        };
    }

    private static JavaPairRDD<IndexPair<String>, Integer> markPairedEnds(final JavaPairRDD<Integer, Iterable<PairedEnds>> keyedPairs,
                                                                          final OpticalDuplicateFinder finder) {
        return keyedPairs.flatMapToPair(keyedPair -> {
            Iterable<PairedEnds> pairedEnds = keyedPair._2();

            final List<Tuple2<IndexPair<String>, Integer>> nonDuplicates = Lists.newArrayList();

            //since we grouped by a non-unique hash code for efficiency we need to regroup by the actual criteria
            //todo this should use library and contig as well probably
            final Collection<List<PairedEnds>> groups = Utils.stream(pairedEnds)
                    .collect(Collectors.groupingBy(pe -> pe.getUnclippedStartPosition() << 8 | pe.getFirstRefIndex())).values();

            for (List<PairedEnds> duplicateGroup : groups) {
                final Map<PairedEnds.Type, List<PairedEnds>> stratifiedByType = splitByType(duplicateGroup);

                // Each key corresponds to either fragments or paired ends, not a mixture of both.
                final List<PairedEnds> fragments = stratifiedByType.get(PairedEnds.Type.FRAGMENT);
                final List<PairedEnds> pairs = stratifiedByType.get(PairedEnds.Type.PAIR);
                final List<PairedEnds> pairsMissingSecondRead = stratifiedByType.get(PairedEnds.Type.PAIRED_BUT_MISSING_SECOND_READ);

                if (fragments != null && !fragments.isEmpty()) { // fragments
                    final Tuple2<IndexPair<String>, Integer> bestFragment = handleFragments(fragments);
                    if( bestFragment != null) {
                        nonDuplicates.add(bestFragment);
                    }
                }

                if (pairs != null && !pairs.isEmpty()) {
                    nonDuplicates.add(handlePairs(pairs, finder));
                }

                if ( pairsMissingSecondRead != null && !pairsMissingSecondRead.isEmpty()){
                    nonDuplicates.addAll(handlePairsMissingSecondRead(pairsMissingSecondRead));
                }
            }

            return nonDuplicates.iterator();
        });
    }

    /**
     * split PairedEnds into groups by their type
     */
    private static Map<PairedEnds.Type, List<PairedEnds>> splitByType(List<PairedEnds> duplicateGroup) {
        final EnumMap<PairedEnds.Type, List<PairedEnds>> byType = new EnumMap<>(PairedEnds.Type.class);
        for(PairedEnds pair: duplicateGroup) {
            byType.compute(pair.getType(), (key, value) -> {
                if (value == null) {
                    final ArrayList<PairedEnds> pairedEnds = new ArrayList<>();
                    pairedEnds.add(pair);
                    return pairedEnds;
                } else {
                    value.add(pair);
                    return value;
                }
            });
        }
        return byType;
    }

    private static List<Tuple2<IndexPair<String>,Integer>> handlePairsMissingSecondRead(List<PairedEnds> pairsMissingSecondRead) {
        // As in Picard, unpaired ends left alone.
        return pairsMissingSecondRead.stream()
                .map(pair -> new Tuple2<>(new IndexPair<>(pair.getName(), pair.getPartitionIndex()), -1))
                .collect(Collectors.toList());
    }

    private static Tuple2<IndexPair<String>, Integer> handlePairs(List<PairedEnds> pairs, OpticalDuplicateFinder finder) {
        final PairedEnds bestPair = pairs.stream()
                .max(PAIRED_ENDS_SCORE_COMPARATOR)
                .orElseThrow(() -> new GATKException.ShouldNeverReachHereException("There was no best pair because the stream was empty, but it shouldn't have been empty."));

        // This must happen last, as findOpticalDuplicates mutates the list.
        // Split by orientation and count duplicates in each group separately.
        // TODO could be better
        final Map<Byte, List<PairedEnds>> groupByOrientation = pairs.stream()
                .peek(pair -> finder.addLocationInformation(pair.getName(), pair))//TODO this needs me to handle the name better
                .collect(Collectors.groupingBy(PairedEnds::getOrientationForOpticalDuplicates));
        final int numOpticalDuplicates;
        //todo do we not have to split the reporting of these by orientation?
        if (groupByOrientation.containsKey(ReadEnds.FR) && groupByOrientation.containsKey(ReadEnds.RF)) {
            final List<PairedEnds> peFR = new ArrayList<>(groupByOrientation.get(ReadEnds.FR));
            final List<PairedEnds> peRF = new ArrayList<>(groupByOrientation.get(ReadEnds.RF));
            numOpticalDuplicates = countOpticalDuplicates(finder, peFR) + countOpticalDuplicates(finder, peRF);
        } else {
            numOpticalDuplicates = countOpticalDuplicates(finder, pairs);
        }
        return (new Tuple2<>(new IndexPair<>(bestPair.getName(), bestPair.getPartitionIndex()), numOpticalDuplicates));
    }

    private static int countOpticalDuplicates(OpticalDuplicateFinder finder, List<PairedEnds> scored) {
        final boolean[] opticalDuplicateFlags = finder.findOpticalDuplicates(scored);
        int numOpticalDuplicates = 0;
        for (final boolean b : opticalDuplicateFlags) {
            if (b) {
                numOpticalDuplicates++;
            }
        }
        return numOpticalDuplicates;
    }

    private static Tuple2<IndexPair<String>, Integer> handleFragments(List<PairedEnds> duplicateFragmentGroup) {
        //empty PairedEnds signify that a pair has a mate somewhere else
        // If there are any non-fragment placeholders at this site, mark everything as duplicates, otherwise compute the best score
        final boolean computeScore = duplicateFragmentGroup.stream().noneMatch(PairedEnds::isEmpty);

        if (computeScore) {
            return duplicateFragmentGroup.stream()
                    .max(PAIRED_ENDS_SCORE_COMPARATOR)
                    .map(best -> new Tuple2<>(new IndexPair<>(best.getName(), best.getPartitionIndex()), -1))
                    .orElse(null);
        } else {
            return null;
        }
    }

    static JavaPairRDD<String, DuplicationMetrics> generateMetrics(final SAMFileHeader header, final JavaRDD<GATKRead> reads) {
        return reads.filter(read -> !read.isSecondaryAlignment() && !read.isSupplementaryAlignment())
                .mapToPair(read -> {
                    final String library = LibraryIdGenerator.getLibraryName(header, read.getReadGroup());
                    DuplicationMetrics metrics = new DuplicationMetrics();
                    metrics.LIBRARY = library;
                    if (read.isUnmapped()) {
                        ++metrics.UNMAPPED_READS;
                    } else if (!read.isPaired() || read.mateIsUnmapped()) {
                        ++metrics.UNPAIRED_READS_EXAMINED;
                    } else {
                        ++metrics.READ_PAIRS_EXAMINED;
                    }

                    if (read.isDuplicate()) {
                        if (!read.isPaired() || read.mateIsUnmapped()) {
                            ++metrics.UNPAIRED_READ_DUPLICATES;
                        } else {
                            ++metrics.READ_PAIR_DUPLICATES;
                        }
                    }
                    if (((SAMRecordToGATKReadAdapter) read).getTransientAttribute(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME)!=null) {
                        metrics.READ_PAIR_OPTICAL_DUPLICATES +=
                                (int)((SAMRecordToGATKReadAdapter) read).getTransientAttribute(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME);
                    }
                    return new Tuple2<>(library, metrics);
                })
                .foldByKey(new DuplicationMetrics(), (metricsSum, m) -> {
                    if (metricsSum.LIBRARY == null) {
                        metricsSum.LIBRARY = m.LIBRARY;
                    }
                    // This should never happen, as we grouped by key using library as the key.
                    if (!metricsSum.LIBRARY.equals(m.LIBRARY)) {
                        throw new GATKException("Two different libraries encountered while summing metrics: " + metricsSum.LIBRARY
                                + " and " + m.LIBRARY);
                    }
                    metricsSum.UNMAPPED_READS += m.UNMAPPED_READS;
                    metricsSum.UNPAIRED_READS_EXAMINED += m.UNPAIRED_READS_EXAMINED;
                    metricsSum.READ_PAIRS_EXAMINED += m.READ_PAIRS_EXAMINED;
                    metricsSum.UNPAIRED_READ_DUPLICATES += m.UNPAIRED_READ_DUPLICATES;
                    metricsSum.READ_PAIR_DUPLICATES += m.READ_PAIR_DUPLICATES;
                    metricsSum.READ_PAIR_OPTICAL_DUPLICATES += m.READ_PAIR_OPTICAL_DUPLICATES;
                    return metricsSum;
                })
                .mapValues(metrics -> {
                    DuplicationMetrics copy = metrics.copy();
                    // Divide these by 2 because they are counted for each read
                    // when they should be counted by pair.
                    copy.READ_PAIRS_EXAMINED = metrics.READ_PAIRS_EXAMINED / 2;
                    copy.READ_PAIR_DUPLICATES = metrics.READ_PAIR_DUPLICATES / 2;

                    copy.calculateDerivedMetrics();
                    if (copy.ESTIMATED_LIBRARY_SIZE == null) {
                        copy.ESTIMATED_LIBRARY_SIZE = 0L;
                    }
                    return copy;
                });
    }

    /**
     * Saves the metrics to a file.
     * Note: the SamFileHeader is needed in order to include libraries that didn't have any duplicates.
     * @param result metrics object, potentially pre-initialized with headers,
     */
    public static void saveMetricsRDD(final MetricsFile<DuplicationMetrics, Double> result, final SAMFileHeader header, final JavaPairRDD<String, DuplicationMetrics> metricsRDD, final String metricsOutputPath) {
        final LibraryIdGenerator libraryIdGenerator = new LibraryIdGenerator(header);

        final Map<String, DuplicationMetrics> nonEmptyMetricsByLibrary = metricsRDD.collectAsMap();           //Unknown Library
        final Map<String, DuplicationMetrics> emptyMapByLibrary = libraryIdGenerator.getMetricsByLibraryMap();//with null

        final List<String> sortedListOfLibraryNames = new ArrayList<>(Sets.union(emptyMapByLibrary.keySet(), nonEmptyMetricsByLibrary.keySet()));
        sortedListOfLibraryNames.sort(Utils.COMPARE_STRINGS_NULLS_FIRST);
        for (final String library : sortedListOfLibraryNames){
            //if a non-empty exists, take it, otherwise take from the the empties. This is done to include libraries with zero data in them.
            //But not all libraries are listed in the header (esp in testing data) so we union empty and non-empty
            final DuplicationMetrics metricsToAdd = nonEmptyMetricsByLibrary.containsKey(library) ? nonEmptyMetricsByLibrary.get(library) : emptyMapByLibrary.get(library);
            metricsToAdd.calculateDerivedMetrics();
            result.addMetric(metricsToAdd);
        }

        if (nonEmptyMetricsByLibrary.size() == 1) {
            result.setHistogram(nonEmptyMetricsByLibrary.values().iterator().next().calculateRoiHistogram());
        }

        MetricsUtils.saveMetrics(result, metricsOutputPath);
    }

    /**
     * GATKRead comparator that compares based on mapping position followed by SAM flags.
     */
    final static class GATKOrder implements Comparator<IndexPair<GATKRead>>, Serializable {
        private static final long serialVersionUID = 1l;
        private final SAMFileHeader header;
        // TODO: Unify with other comparators in the codebase

        public GATKOrder(final SAMFileHeader header) {
            this.header = header;
        }

        @Override
        public int compare(final IndexPair<GATKRead> lhsPair, final IndexPair<GATKRead> rhsPair) {
            final GATKRead lhs = lhsPair.getValue();
            final GATKRead rhs = rhsPair.getValue();
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

            //Note: negate the result because we want first-of-pair to be before second
            //ie, want 'second' to be sorted after first, so want to return -1 for (true, false)
            final int res7 = -Boolean.compare(lhs.isFirstOfPair(), rhs.isFirstOfPair());
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

    /**
     * Comparator for sorting Reads by coordinate. Note that a header is required in
     * order to meaningfully compare contigs.
     *
     * Uses the various other fields in a read to break ties for reads that share
     * the same location.
     *
     * Ordering is almost identical to the {@link htsjdk.samtools.SAMRecordCoordinateComparator},
     * modulo a few subtle differences in tie-breaking rules for reads that share the same
     * position. This comparator will produce an ordering consistent with coordinate ordering
     * in a bam file, including interleaving unmapped reads assigned the positions of their
     * mates with the mapped reads.
     */
    public static final class PairedEndsCoordinateComparator implements Comparator<PairedEnds>, Serializable {
        private static final long serialVersionUID = 1L;

        public static final PairedEndsCoordinateComparator INSTANCE = new PairedEndsCoordinateComparator();
        private PairedEndsCoordinateComparator() { }

        @Override
        public int compare( PairedEnds first, PairedEnds second ) {
            int result = compareCoordinates(first, second);
            if ( result != 0 ) {
                return result;
            }

            //This is done to mimic SAMRecordCoordinateComparator's behavior
            if (first.isR1R() != second.isR1R()) {
                return first.isR1R() ? -1: 1;
            }

            if ( first.getName() != null && second.getName() != null ) {
                result = first.getName().compareTo(second.getName());
            }
            return result;
        }

        public static int compareCoordinates(final PairedEnds first, final PairedEnds second ) {
            final int firstRefIndex = first.getFirstRefIndex();
            final int secondRefIndex = second.getFirstRefIndex();

            if ( firstRefIndex == -1 ) {
                return (secondRefIndex == -1 ? 0 : 1);
            }
            else if ( secondRefIndex == -1 ) {
                return -1;
            }

            final int refIndexDifference = firstRefIndex - secondRefIndex;
            if ( refIndexDifference != 0 ) {
                return refIndexDifference;
            }

            return Integer.compare(first.getFirstStartPosition(), second.getFirstStartPosition());
        }
    }
}