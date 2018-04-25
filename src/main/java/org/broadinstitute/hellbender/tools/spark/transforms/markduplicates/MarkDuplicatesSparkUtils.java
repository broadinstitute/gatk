package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.serializers.FieldSerializer;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.markduplicates.*;
import org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords.*;
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
    // This comparator represents the tiebreaking for PairedEnds duplicate marking.
    // We compare first on score, followed by unclipped start position (which is reversed here because of the expected ordering)
    private static final Comparator<PairedEnds> PAIRED_ENDS_SCORE_COMPARATOR = Comparator.comparing(PairedEnds::getScore)
            .thenComparing(PairedEndsCoordinateComparator.INSTANCE.reversed());

    /**
     * Wrapper object used for storing an object and some type of index information.
     *
     * MarkDuplicates uses this object remember which partition of the original bam each read came from in order to
     * efficiently zip their duplicate marked data back into the correct place without shuffling the original bam.
     */
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

        @VisibleForTesting
        IndexPair(T value, int index) {
            this.value = value;
            this.index = index;
        }
    }

    /**
     * (0) filter: remove unpaired reads and reads with an unmapped mate.
     * (1) keyReadsByName: label each read with its read group and read name.
     * (2) GroupByKey: group together reads with the same group and name.
     * (3) keyMarkDuplicatesSparkRecords with alignment info:
     *   (a) Generate a fragment or emptyFragment from each read if it's unpaired.
     *   (b) Pair grouped reads into MarkDuplicatesSparkRecord. In most cases there will only be two reads
     *       with the same name. Mapped reads missing mates will be emitted as fragments, more than two reads will cause an exception.
     *   (c) Label each read with alignment information: Library, reference index,
     *       stranded unclipped start and reverse strand.
     *   (d) Unmapped Pairs, Templates of entirely non-primary reads, etc are passed through as unmarked reads
     * (4) GroupByKey: Group MarkDuplicatesSparkRecord that share alignment information. These pairs
     *     are duplicates of each other.
     * (5) markDuplicatePairs:
     *   (a) For each group created by (4), sort the pairs by score and mark all but the
     *       highest scoring as duplicates.
     *   (b) Determine which duplicates are optical duplicates and increase the overall count.
     */
    static JavaPairRDD<IndexPair<String>, Integer> transformToDuplicateNames(final SAMFileHeader header, final MarkDuplicatesScoringStrategy scoringStrategy, final OpticalDuplicateFinder finder, final JavaRDD<GATKRead>  reads, final int numReducers) {
        // we treat these specially and don't mark them as duplicates
        final JavaRDD<GATKRead> mappedReads = reads.filter(ReadFilterLibrary.MAPPED::test);

        final JavaPairRDD<String, Iterable<IndexPair<GATKRead>>> keyedReads = getReadsGroupedByName(header, mappedReads, numReducers);

        // Place all the reads into a single RDD of MarkDuplicatesSparkRecord objects
        final JavaPairRDD<Integer, MarkDuplicatesSparkRecord> pairedEnds = keyedReads.flatMapToPair(keyedRead -> {
            final List<Tuple2<Integer, MarkDuplicatesSparkRecord>> out = Lists.newArrayList();
            final IndexPair<?>[] hadNonPrimaryRead = {null};

            final List<IndexPair<GATKRead>> primaryReads = Utils.stream(keyedRead._2())
                    ////// Making The Fragments //////
                    // Make a PairedEnd object with no second read for each fragment (and an empty one for each paired read)
                    .peek(readWithIndex -> {
                        final GATKRead read = readWithIndex.getValue();
                        if (!(read.isSecondaryAlignment()||read.isSupplementaryAlignment())) {
                            PairedEnds fragment = (ReadUtils.readHasMappedMate(read)) ?
                                    MarkDuplicatesSparkRecord.newEmptyFragment(read, header) :
                                    MarkDuplicatesSparkRecord.newFragment(read, header, readWithIndex.getIndex(), scoringStrategy);

                            out.add(new Tuple2<>(fragment.key(), fragment));
                        } else {
                            hadNonPrimaryRead[0] = readWithIndex;
                        }
                    })
                    .filter(readWithIndex -> ReadUtils.readHasMappedMate(readWithIndex.getValue()))

            ////// Making The Paired Reads //////
            // Write each paired read with a mapped mate as a pair
                    .filter(indexPair -> !(indexPair.getValue().isSecondaryAlignment()||indexPair.getValue().isSupplementaryAlignment()))
                    .collect(Collectors.toList());

            // Mark duplicates cant properly handle templates with more than two reads in a pair
            if (primaryReads.size()>2) {
                throw new UserException.UnimplementedFeature(String.format("MarkDuplicatesSpark only supports singleton fragments and pairs. We found the following group with >2 primary reads: ( %d number of reads)." +
                        " \n%s.", primaryReads.size(),primaryReads.stream().map(Object::toString).collect(Collectors.joining("\n"))));

            // If there are two primary reads in the group pass them as a pair
            } else if (primaryReads.size()==2) {
                final IndexPair<GATKRead> firstRead = primaryReads.get(0);
                final IndexPair<GATKRead> secondRead = primaryReads.get(1);
                final PairedEnds pair = MarkDuplicatesSparkRecord.newPair(firstRead.getValue(), secondRead.getValue(), header, secondRead.getIndex(), scoringStrategy);
                out.add(new Tuple2<>(pair.key(), pair));

            // If there is one paired read in the template this probably means the bam is missing its mate, don't duplicate mark it
            } else if (primaryReads.size()==1) {
                final IndexPair<GATKRead> firstRead = primaryReads.get(0);
                final MarkDuplicatesSparkRecord pass = MarkDuplicatesSparkRecord.getPassthrough(firstRead.getValue(), firstRead.getIndex());
                out.add(new Tuple2<>(pass.key(), pass));

            // else that means there are no non-secondary or supplementary reads, thus we want the group to pass through through unmarked
            } else {
                if (hadNonPrimaryRead[0] !=null) {
                    final MarkDuplicatesSparkRecord pass = MarkDuplicatesSparkRecord.getPassthrough((GATKRead)hadNonPrimaryRead[0].getValue(), hadNonPrimaryRead[0].getIndex());
                    out.add(new Tuple2<>(pass.key(), pass));
                }
            }

            return out.iterator();
        });

        final JavaPairRDD<Integer, Iterable<MarkDuplicatesSparkRecord>> keyedPairs = pairedEnds.groupByKey(); //TODO evaluate replacing this with a smart aggregate by key.

        return markDuplicateRecords(keyedPairs, finder);
    }

    /**
     * Method that ensures the reads are grouped together keyed by their readname groups with indexpairs represeting the source partition.
     * If the bam is querygrouped/queryname sorted then it calls spanReadsByKey to perform the mapping operation
     * If the bam is sorted in some other way it performs a groupBy operation on the key
     */
    private static JavaPairRDD<String, Iterable<IndexPair<GATKRead>>> getReadsGroupedByName(SAMFileHeader header, JavaRDD<GATKRead> reads, int numReducers) {

        final JavaPairRDD<String, Iterable<IndexPair<GATKRead>>> keyedReads;
        final JavaRDD<IndexPair<GATKRead>> indexedReads = reads.mapPartitionsWithIndex(
                (index, iter) -> Utils.stream(iter).map(read -> {
                    if (!(read.getClass() == SAMRecordToGATKReadAdapter.class)) {
                        throw new GATKException(String.format("MarkDuplicatesSpark currently only supports SAMRecords as an underlying reads data source class, %s found instead",
                                read.getClass().toString()));
                }
                return new IndexPair<>(read, index);}).iterator(), false);
        if (ReadUtils.isReadNameGroupedBam(header)) {
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

    /**
     * Method which takes an RDD of reads that is guaranteed to have every readname group placed together on the same
     * partition and maps those so a JavaPairRDD with the readname as the key.
     */
    private static JavaPairRDD<String, Iterable<IndexPair<GATKRead>>> spanReadsByKey(final JavaRDD<IndexPair<GATKRead>> reads) {
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

    /**
     * Primary landing point for MarkDuplicateSparkRecords:
     *  - Handles separating out hashed keys into into groups by start position/readgroup
     *  - Further separates out MarkDuplicatesSparkRecord by their record objects
     *  - Farms out to methods which handles each group
     *  - Collects the results and returns an iterator
     */
    @SuppressWarnings("unchecked")
    private static JavaPairRDD<IndexPair<String>, Integer> markDuplicateRecords(final JavaPairRDD<Integer, Iterable<MarkDuplicatesSparkRecord>> keyedPairs,
                                                                                final OpticalDuplicateFinder finder) {
        return keyedPairs.flatMapToPair(keyedPair -> {
            Iterable<MarkDuplicatesSparkRecord> pairGroups = keyedPair._2();

            final List<Tuple2<IndexPair<String>, Integer>> nonDuplicates = Lists.newArrayList();

            //since we grouped by a non-unique hash code for efficiency we need to regroup by the actual criteria
            //todo this should use library and contig as well probably
            //todo these should all be one traversal over the records)
            final Collection<List<MarkDuplicatesSparkRecord>> groups = Utils.stream(pairGroups)
                    .collect(Collectors.groupingBy(MarkDuplicatesSparkUtils::getGroupKey)).values();

            for (List<MarkDuplicatesSparkRecord> duplicateGroup : groups) {
                final Map<MarkDuplicatesSparkRecord.Type, List<MarkDuplicatesSparkRecord>> stratifiedByType = splitByType(duplicateGroup);

                // Each key corresponds to either fragments or paired ends, not a mixture of both.
                final List<MarkDuplicatesSparkRecord> emptyFragments = stratifiedByType.get(MarkDuplicatesSparkRecord.Type.EMPTY_FRAGMENT);
                final List<MarkDuplicatesSparkRecord> fragments = stratifiedByType.get(MarkDuplicatesSparkRecord.Type.FRAGMENT);
                final List<Pair> pairs = (List<Pair>) (List)(stratifiedByType.get(MarkDuplicatesSparkRecord.Type.PAIR));
                final List<MarkDuplicatesSparkRecord> passthroughs = stratifiedByType.get(MarkDuplicatesSparkRecord.Type.PASSTHROUGH);

                //empty MarkDuplicatesSparkRecord signify that a pair has a mate somewhere else
                // If there are any non-fragment placeholders at this site, mark everything as duplicates, otherwise compute the best score
                if (Utils.isNonEmpty(fragments) && !Utils.isNonEmpty(emptyFragments)) {
                    final Tuple2<IndexPair<String>, Integer> bestFragment = handleFragments(fragments);
                    nonDuplicates.add(bestFragment);
                }

                if (Utils.isNonEmpty(pairs)) {
                    nonDuplicates.add(handlePairs(pairs, finder));
                }

                if (Utils.isNonEmpty(passthroughs)) {
                    nonDuplicates.addAll(handlePassthroughs(passthroughs));
                }
            }

            return nonDuplicates.iterator();
        });
    }

    // Note, this uses bitshift operators in order to perform only a single groupBy operation for all the merged data
    private static long getGroupKey(MarkDuplicatesSparkRecord record) {
        return record.getClass()==Passthrough.class?-1:
                (((long)((PairedEnds)record).getUnclippedStartPosition()) << 32 |
                        ((PairedEnds)record).getFirstRefIndex() << 16 );
        //| ((PairedEnds)pe).getLibraryIndex())).values();
    }

    /**
     * split MarkDuplicatesSparkRecord into groups by their type
     */
    private static Map<MarkDuplicatesSparkRecord.Type, List<MarkDuplicatesSparkRecord>> splitByType(List<MarkDuplicatesSparkRecord> duplicateGroup) {
        final EnumMap<MarkDuplicatesSparkRecord.Type, List<MarkDuplicatesSparkRecord>> byType = new EnumMap<>(MarkDuplicatesSparkRecord.Type.class);
        for(MarkDuplicatesSparkRecord pair: duplicateGroup) {
            byType.compute(pair.getType(), (key, value) -> {
                if (value == null) {
                    final ArrayList<MarkDuplicatesSparkRecord> pairedEnds = new ArrayList<>();
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

        private static List<Tuple2<IndexPair<String>,Integer>> handlePassthroughs(List<MarkDuplicatesSparkRecord> passthroughs) {
        // Emit the passthrough reads as non-duplicates.
        return passthroughs.stream()
                .map(pair -> new Tuple2<>(new IndexPair<>(pair.getName(), pair.getPartitionIndex()), -1))
                .collect(Collectors.toList());
    }

    private static Tuple2<IndexPair<String>, Integer> handlePairs(List<Pair> pairs, OpticalDuplicateFinder finder) {
        final MarkDuplicatesSparkRecord bestPair = pairs.stream()
                .max(PAIRED_ENDS_SCORE_COMPARATOR)
                .orElseThrow(() -> new GATKException.ShouldNeverReachHereException("There was no best pair because the stream was empty, but it shouldn't have been empty."));

        // Split by orientation and count duplicates in each group separately.
        final Map<Byte, List<Pair>> groupByOrientation = pairs.stream()
                .peek(pair -> finder.addLocationInformation(pair.getName(), pair))//TODO this needs me to handle the name better
                .collect(Collectors.groupingBy(Pair::getOrientationForOpticalDuplicates));
        final int numOpticalDuplicates;
        //todo do we not have to split the reporting of these by orientation?
        if (groupByOrientation.containsKey(ReadEnds.FR) && groupByOrientation.containsKey(ReadEnds.RF)) {
            final List<Pair> peFR = new ArrayList<>(groupByOrientation.get(ReadEnds.FR));
            final List<Pair> peRF = new ArrayList<>(groupByOrientation.get(ReadEnds.RF));
            numOpticalDuplicates = countOpticalDuplicates(finder, peFR) + countOpticalDuplicates(finder, peRF);
        } else {
            numOpticalDuplicates = countOpticalDuplicates(finder, pairs);
        }
        return (new Tuple2<>(new IndexPair<>(bestPair.getName(), bestPair.getPartitionIndex()), numOpticalDuplicates));
    }

    private static int countOpticalDuplicates(OpticalDuplicateFinder finder, List<Pair> scored) {
        final boolean[] opticalDuplicateFlags = finder.findOpticalDuplicates(scored);
        int numOpticalDuplicates = 0;
        for (final boolean b : opticalDuplicateFlags) {
            if (b) {
                numOpticalDuplicates++;
            }
        }
        return numOpticalDuplicates;
    }

    /**
     * If there are fragments with no non-fragments overlapping at a site, select the best one according to PAIRED_ENDS_SCORE_COMPARATOR
     */
    private static Tuple2<IndexPair<String>, Integer> handleFragments(List<MarkDuplicatesSparkRecord> duplicateFragmentGroup) {
        return duplicateFragmentGroup.stream()
                    .map(f -> (Fragment)f)
                    .max(PAIRED_ENDS_SCORE_COMPARATOR)
                    .map(best -> new Tuple2<>(new IndexPair<>(best.getName(), best.getPartitionIndex()), -1))
                    .orElse(null);
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
                    // NOTE: we use the SAMRecord transientAttribute field here specifically to prevent the already
                    // serialized read from being parsed again here for performance reasons.
                    if (((SAMRecordToGATKReadAdapter) read).getTransientAttribute(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME)!=null) {
                        // NOTE: there is a safety check above in getReadsGroupedByName()
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
        for (final String library : sortedListOfLibraryNames) {
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
     * Comparator for sorting Reads by coordinate. Note that a header is required in
     * order to meaningfully compare contigs.
     *
     * Uses the various other fields in a read to break ties for reads that share
     * the same location.
     *
     * Ordering is not the same as {@link htsjdk.samtools.SAMRecordCoordinateComparator}.
     * It compares two pairedEnds objects by their first reads according to first their clipped start positions
     * (they were matched together based on UnclippedStartOriginally), then the orientation of the strand, followed by
     * the readname lexicographical order.
     *
     * NOTE: Because the original records were grouped by readname, we know that they must be unique once they hit this
     *       comparator and thus we don't need to worry about further tiebreaking for this method.
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