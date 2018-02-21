package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import com.google.common.collect.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.*;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Utility classes and functions for Mark Duplicates.
 */
public class MarkDuplicatesSparkUtils {

    // Used to set an attribute on the GATKRead marking this read as an optical duplicate.
    public static final String OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME = "OD";

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
    static JavaPairRDD<String, Integer> transformToDuplicateNames(final SAMFileHeader header, final MarkDuplicatesScoringStrategy scoringStrategy, final OpticalDuplicateFinder finder, final JavaRDD<GATKRead>  reads, final int numReducers) {

        JavaPairRDD<String, Iterable<GATKRead>> keyedReads;
        if (SAMFileHeader.SortOrder.queryname.equals(header.getSortOrder())) {
            // reads are already sorted by name, so perform grouping within the partition (no shuffle)
            keyedReads = spanReadsByKey(header, reads);
        } else {
            // sort by group and name (incurs a shuffle)
            JavaPairRDD<String, GATKRead> keyReadPairs = reads.mapToPair(read -> new Tuple2<>(ReadsKey.keyForRead(header, read), read));
            keyedReads = keyReadPairs.groupByKey(numReducers);
        }

        // Place all the reads into a single RDD seperated
        JavaPairRDD<Integer, Iterable<PairedEnds>> keyedPairs = keyedReads.flatMapToPair(keyedRead -> {
            final List<Tuple2<Integer, PairedEnds>> out = Lists.newArrayList();

            final Iterator<GATKRead> sorted = StreamSupport.stream(keyedRead._2().spliterator(), false)
                    ////// Making The Fragments //////
                    // Make a PairedEnd object with no second read for each fragment (and an empty one for each paired read)
                    .peek(read -> {
                        //read.setIsDuplicate(false);

                        PairedEnds fragment = (ReadUtils.readHasMappedMate(read)) ?
                                PairedEnds.empty(ReadUtils.getStrandedUnclippedStart(read)) :
                                PairedEnds.of(read, header);

                        out.add(new Tuple2<>(fragment.isEmpty()?
                                            ReadsKey.keyForFragment(header, read, fragment.getStartPosition()) :
                                            fragment.keyForFragment(header),
                                fragment));
                        if (!fragment.isEmpty()) fragment.score(scoringStrategy);
                    })
                    .filter(ReadUtils::readHasMappedMate)
                    .sorted(new GATKOrder(header))
                    .iterator();

            ////// Making The Paired Reads //////
            // Write each paired read with a mapped mate as a pair
            PairedEnds pair;
            while (sorted.hasNext()) {
                pair = PairedEnds.of(sorted.next(), false, header);
                if (sorted.hasNext()) {
                    pair.and(sorted.next(), header);
                    out.add(new Tuple2<>(pair.key(header), pair));
                    pair.score(scoringStrategy);
                }
                else {
                    out.add(new Tuple2<>(pair.key(header), pair));
                    pair.score(scoringStrategy);
                }
            }
            return out.iterator();
        }).groupByKey(numReducers);

        return markPairedEnds(keyedPairs, scoringStrategy, finder, header);
    }

    static JavaPairRDD<String, Iterable<GATKRead>> spanReadsByKey(final SAMFileHeader header, final JavaRDD<GATKRead> reads) {
        JavaPairRDD<String, GATKRead> nameReadPairs = reads.mapToPair(read -> new Tuple2<>(read.getName(), read));
        return spanByKey(nameReadPairs).flatMapToPair(namedRead -> {
            // for each name, separate reads by key (group name)
            List<Tuple2<String, Iterable<GATKRead>>> out = Lists.newArrayList();
            ListMultimap<String, GATKRead> multi = LinkedListMultimap.create();
            for (GATKRead read : namedRead._2()) {
                multi.put(ReadsKey.keyForRead(header, read), read);
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
    static <K, V> JavaPairRDD<K, Iterable<V>> spanByKey(JavaPairRDD<K, V> rdd) {
        return rdd.mapPartitionsToPair(iter -> spanningIterator(iter));
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

//    private static class ReadMinimalRepresentation {
//        String name;
//        int readIndex;
//        int startPosition;
//        int unclippedStartPostion;
//
//        ReadMinimalRepresentation(final GATKRead read, final SAMFileHeader header) {
//            this.name = read.getName();
//            this.readIndex = ReadUtils.getReferenceIndex(read, header);
//        }
//
//        static ReadMinimalRepresentation getRepresntation(final GATKRead read, final SAMFileHeader header) {
//            return new ReadMinimalRepresentation(read, header);
//        }
//    }


    static JavaPairRDD<String, Integer> markPairedEnds(final JavaPairRDD<Integer, Iterable<PairedEnds>> keyedPairs,
                                            final MarkDuplicatesScoringStrategy scoringStrategy,
                                            final OpticalDuplicateFinder finder, final SAMFileHeader header) {
        return keyedPairs.flatMapToPair(keyedPair -> {
            Iterable<PairedEnds> pairedEnds = keyedPair._2();
            final ImmutableListMultimap<Boolean, PairedEnds> stratifiedByFragments = Multimaps.index(pairedEnds, PairedEnds::isFragment);

            Comparator<PairedEnds> pairedEndsComparator =
                    Comparator.<PairedEnds, Integer>comparing(pe -> pe.score(scoringStrategy)).reversed()
                            .thenComparing((o1, o2) -> new ReadCoordinateComparator(header).compare(o1.first(), o2.first()));

            // Each key corresponds to either fragments or paired ends, not a mixture of both.

            final List<Tuple2<String, Integer>> out = Lists.newArrayList();

            /////// fragments ////////
            if (!stratifiedByFragments.get(true).isEmpty()) { // fragments
                Tuple2<String, Integer> fragmentNonDuplicate = handleFragments(stratifiedByFragments.get(true), pairedEndsComparator);
                if (fragmentNonDuplicate!=null) {
                    out.add(fragmentNonDuplicate);
                }
            }

            /////// ReadPairs ////////
            if (!stratifiedByFragments.get(false).isEmpty()) {
                final ImmutableListMultimap<Integer, PairedEnds> groups = Multimaps.index(stratifiedByFragments.get(false), PairedEnds::getStartPosition);

                for (int key : groups.keys()) {
//                    // As in Picard, unpaired ends left alone.
//                    groups.get(key).stream().filter(pair -> pair.second() == null).forEach(pair -> out.add(pair.first()));
//
//                    // Order by score using ReadCoordinateComparator for tie-breaking.
//                    final List<PairedEnds> scored = groups.get(key).stream().filter(pair -> pair.second() != null).sorted(pairedEndsComparator).collect(Collectors.toList());
//
//                    final PairedEnds best = Iterables.getFirst(scored, null);
//                    if (best == null) {
//                        return out.iterator();
//                    }

                    PairedEnds maxScore = null;
                    for (PairedEnds pair : groups.get(key)) {
                        GATKRead read = pair.first();

                        // As in Picard, unpaired ends left alone.
                        if (pair.second() == null) {
                            out.add(new Tuple2<String,Integer>(pair.getName(),-1));

                            // Order by score using ReadCoordinateComparator for tie-breaking.
                        } else {

//                            out.add(pair.first());
//                            out.add(pair.second());

                            // There are no paired reads, mark all but the highest scoring fragment as duplicate.
                            if (maxScore != null) {
                                maxScore = (pairedEndsComparator.compare(pair, maxScore) > 0) ? maxScore : pair;
                            } else {
                                maxScore = pair;
                            }

                            finder.addLocationInformation(pair.first().getName(), pair);//TODO this needs me to handle the name better
                        }
                    }
                    if (maxScore == null) {
                        return out.iterator();
                    }
//
//                    // Mark everyone who's not best as a duplicate
//                    for (final PairedEnds pair : Iterables.skip(scored, 1)) {
//                        pair.first().setIsDuplicate(true);
//                        pair.second().setIsDuplicate(true);
//                    }
//
//                    // Now, add location information to the paired ends
//                    for (final PairedEnds pair : scored) {
//                        // Both elements in the pair have the same name
//                        finder.addLocationInformation(pair.first().getName(), pair);
//                    }

                    // This must happen last, as findOpticalDuplicates mutates the list.
                    // Split by orientation and count duplicates in each group separately.
                    // TODO could be better
                    List<PairedEnds> scored = stratifiedByFragments.get(false).stream().filter(pair -> pair.second() != null).collect(Collectors.toList());
                    final ImmutableListMultimap<Byte, PairedEnds> groupByOrientation = Multimaps.index(scored, pe -> pe.getOrientationForOpticalDuplicates());
                    final int numOpticalDuplicates;
                    if (groupByOrientation.containsKey(ReadEnds.FR) && groupByOrientation.containsKey(ReadEnds.RF)) {
                        final List<PairedEnds> peFR = new ArrayList<>(groupByOrientation.get(ReadEnds.FR));
                        final List<PairedEnds> peRF = new ArrayList<>(groupByOrientation.get(ReadEnds.RF));
                        numOpticalDuplicates = countOpticalDuplicates(finder, peFR) + countOpticalDuplicates(finder, peRF);
                    } else {
                        numOpticalDuplicates = countOpticalDuplicates(finder, scored);
                    }
                    //maxScore.first().setAttribute(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME, numOpticalDuplicates);
                    out.add(new Tuple2<>(maxScore.getName(), numOpticalDuplicates));
                }
            }
            return out.iterator();
        });
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

    private static Tuple2<String,Integer> handleFragments(Iterable<PairedEnds> pairedEnds, Comparator<PairedEnds> comparator) {

        final ImmutableListMultimap<Integer, PairedEnds> groups = Multimaps.index(pairedEnds, pair -> pair.getStartPosition());

        for (int key : groups.keys()) {
            final Map<Boolean, List<PairedEnds>> byPairing = groups.get(key).stream()
                    .collect(Collectors.partitioningBy(
                            PairedEnds::isEmpty
                    ));

            // If there are any non-fragments at this site, mark everything as duplicates, otherwise compute the best score
            //TODO this must be collapsed 
            boolean computeScore = byPairing.get(true).isEmpty();
            PairedEnds maxScore = null;
            for (PairedEnds pairing : byPairing.get(false)) {

                // There are no paired reads, mark all but the highest scoring fragment as duplicate.
                if (computeScore) {
                    if (maxScore != null) {
                        maxScore = (comparator.compare(pairing, maxScore) > 0) ? maxScore : pairing;
                    } else {
                        maxScore = pairing;
                    }
                }

            }

            // Mark best scoring read as a non-duplicate
            if (maxScore != null) {
               return new Tuple2<String,Integer>(maxScore.getName(),-1);
            }
        }
        return null;

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
                    if (read.hasAttribute(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME)) {
                        metrics.READ_PAIR_OPTICAL_DUPLICATES +=
                                read.getAttributeAsInteger(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME);
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
}