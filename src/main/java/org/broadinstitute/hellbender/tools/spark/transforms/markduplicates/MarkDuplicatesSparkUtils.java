package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.serializers.FieldSerializer;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapGroupsFunction;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.sql.*;
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
import org.broadinstitute.hellbender.utils.read.markduplicates.sparksql.NonDuplicatesSparkSqlRecord;
import org.broadinstitute.hellbender.utils.read.markduplicates.sparksql.MarkDuplicatesSparkSqlRecord;
import org.broadinstitute.hellbender.utils.read.markduplicates.sparksql.ReadsKeySparkSqlRecord;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import picard.sam.markduplicates.util.ReadEnds;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Utility classes and functions for Mark Duplicates.
 */
public class MarkDuplicatesSparkUtils {

    // Used to set an attribute on the GATKRead marking this read as an optical duplicate.
    public static final String OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME = "OD";
    // This comparator represents the tiebreaking for PairedEnds duplicate marking.
    // We compare first on score, followed by unclipped start position (which is reversed here because of the expected ordering)
    private static final Comparator<MarkDuplicatesSparkSqlRecord> PAIRED_ENDS_SCORE_COMPARATOR = Comparator.comparing(MarkDuplicatesSparkSqlRecord::getScore)
            .thenComparing(MarkDuplicatesSparkSqlRecordPhysicalLocationComparator.INSTANCE.reversed());

    /**
     * Returns the library associated with the provided read's read group.
     * Or the specified default if no library is found
     *
     * @param read read whose library to retrieve
     * @param header SAM header containing read groups
     * @return the library for the provided read's read group as a String,
     *         or the default value if the read has no read group.
     */
    public static String getLibraryForRead(final GATKRead read, final SAMFileHeader header, String defaultLibrary) {
        final SAMReadGroupRecord readGroup = ReadUtils.getSAMReadGroupRecord(read, header);
        if (readGroup != null) {
            String library = readGroup.getLibrary();
            return library==null? defaultLibrary : library;
        } else {
            if (read.getReadGroup() == null) {
                throw new UserException.ReadMissingReadGroup(read);
            } else {
                throw new UserException.HeaderMissingReadGroup(read);
            }
        }
    }

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

        @Override
        public String toString() {
            return "indexpair["+index+","+value.toString()+"]";
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
    static JavaPairRDD<IndexPair<String>, Integer> transformToDuplicateNames(final SAMFileHeader header, final MarkDuplicatesScoringStrategy scoringStrategy, final OpticalDuplicateFinder finder, final JavaRDD<GATKRead>  reads, final int numReducers, final boolean markOpticalDups) {
        // we treat these specially and don't mark them as duplicates
        final JavaRDD<GATKRead> mappedReads = reads.filter(ReadFilterLibrary.MAPPED::test);

        final JavaPairRDD<String, Iterable<IndexPair<GATKRead>>> keyedReads = getReadsGroupedByName(header, mappedReads, numReducers);

        final Broadcast<Map<String, Short>> headerReadGroupIndexMap = JavaSparkContext.fromSparkContext(reads.context()).broadcast( getHeaderReadGroupIndexMap(header));
        final Broadcast<Map<String, Byte>> libraryIndex = JavaSparkContext.fromSparkContext(reads.context()).broadcast( constructLibraryIndex(header));

        // Place all the reads into a single RDD of MarkDuplicatesSparkSqlRecord objects
        final JavaRDD<MarkDuplicatesSparkSqlRecord> pairedEnds = keyedReads.flatMap(keyedRead -> {
            final List<MarkDuplicatesSparkSqlRecord> out = Lists.newArrayList();
            final IndexPair<?>[] hadNonPrimaryRead = {null};

            final List<IndexPair<GATKRead>> primaryReads = Utils.stream(keyedRead._2())
                    ////// Making The Fragments //////
                    // Make a PairedEnd object with no second read for each fragment (and an empty one for each paired read)
                    .peek(readWithIndex -> {
                        final GATKRead read = readWithIndex.getValue();
                        if (!(read.isSecondaryAlignment()||read.isSupplementaryAlignment())) {
                            MarkDuplicatesSparkSqlRecord fragment = (ReadUtils.readHasMappedMate(read)) ?
                                    MarkDuplicatesSparkSqlRecord.newEmptyFragment(read, header, libraryIndex.getValue()) :
                                    MarkDuplicatesSparkSqlRecord.newFragment(read, header, readWithIndex.getIndex(), scoringStrategy, libraryIndex.getValue());

                            out.add(fragment);
                        } else {
                            hadNonPrimaryRead[0] = readWithIndex;
                        }
                    })
                    .filter(indexPair -> !(indexPair.getValue().isSecondaryAlignment()||indexPair.getValue().isSupplementaryAlignment()))
                    .collect(Collectors.toList());

            // Catching the case where there are only secondary and supplementary reads in the readname group
            if (primaryReads.isEmpty()) {
                final MarkDuplicatesSparkSqlRecord pass = MarkDuplicatesSparkSqlRecord.getPassthrough((GATKRead)hadNonPrimaryRead[0].getValue(), hadNonPrimaryRead[0].getIndex());
                out.add(pass);
                return out.iterator();

                // Mark duplicates cant properly handle templates with more than two reads in a pair
            } else if (primaryReads.size()>2) {
                throw new UserException.UnimplementedFeature(String.format("MarkDuplicatesSpark only supports singleton fragments and pairs. We found the following group with >2 primary reads: ( %d number of reads)." +
                        " \n%s.", primaryReads.size(), primaryReads.stream().map(Object::toString).collect(Collectors.joining("\n"))));
            }

            ////// Making The Paired Reads //////
            // Write each paired read with a mapped mate as a pair
            final List<IndexPair<GATKRead>> mappedPair = primaryReads.stream()
                    .filter(readWithIndex -> ReadUtils.readHasMappedMate(readWithIndex.getValue()))
                    .collect(Collectors.toList());

            // If there are two primary reads in the group pass them as a pair
            if (mappedPair.size()==2) {
                final GATKRead firstRead = mappedPair.get(0).getValue();
                final IndexPair<GATKRead> secondRead = mappedPair.get(1);
                final MarkDuplicatesSparkSqlRecord pair = MarkDuplicatesSparkSqlRecord.newPair(firstRead, secondRead.getValue(), header, secondRead.getIndex(), scoringStrategy, libraryIndex.getValue());
                // Validate and add the read group to the pair
                final Short readGroup = headerReadGroupIndexMap.getValue().get(firstRead.getReadGroup());
                if (readGroup != null) {
                    pair.setReadGroup(readGroup);
                } else {
                    throw (firstRead.getReadGroup()==null) ?
                            new UserException.ReadMissingReadGroup(firstRead) :
                            new UserException.HeaderMissingReadGroup(firstRead);
                }
                out.add(pair);

                // If there is one paired read in the template this probably means the bam is missing its mate, don't duplicate mark it
            } else if (mappedPair.size()==1) {
                final IndexPair<GATKRead> firstRead = mappedPair.get(0);
                final MarkDuplicatesSparkSqlRecord pass = MarkDuplicatesSparkSqlRecord.getPassthrough(firstRead.getValue(), firstRead.getIndex());
                out.add(pass);
            }
            // If mappedPair is empty here, it probably means that we had a fragment with an unmapped mate, which has already been built
            // and added to out. So we just pass through and return.

            return out.iterator();
        });

        // convert to Spark Dataset
        SparkSession sparkSession = SparkSession.builder().sparkContext(JavaSparkContext.fromSparkContext(reads.context()).sc()).getOrCreate();
        Dataset<MarkDuplicatesSparkSqlRecord> markDuplicatesSparkSqlRecordDataset = sparkSession.createDataset(pairedEnds.rdd(), Encoders.bean(MarkDuplicatesSparkSqlRecord.class));
        // group by key
        KeyValueGroupedDataset<ReadsKeySparkSqlRecord, MarkDuplicatesSparkSqlRecord> grouped = markDuplicatesSparkSqlRecordDataset.groupByKey(
                (MapFunction<MarkDuplicatesSparkSqlRecord, ReadsKeySparkSqlRecord>) MarkDuplicatesSparkSqlRecord::key,
                Encoders.bean(ReadsKeySparkSqlRecord.class));
        // and mark duplicates
        Dataset<NonDuplicatesSparkSqlRecord> nonDuplicates = grouped.flatMapGroups(
                (FlatMapGroupsFunction<ReadsKeySparkSqlRecord, MarkDuplicatesSparkSqlRecord, NonDuplicatesSparkSqlRecord>) (readsKeySparkSqlRecord, iterator) ->
                        markDuplicateRecords(Utils.stream(iterator).collect(Collectors.toList()), finder, markOpticalDups), Encoders.bean(NonDuplicatesSparkSqlRecord.class));
        // convert back to RDD
        return nonDuplicates.toJavaRDD().mapToPair((PairFunction<NonDuplicatesSparkSqlRecord, IndexPair<String>, Integer>) record ->
                new Tuple2<>(new IndexPair<>(record.getName(), record.getPartitionIndex()), record.getNumOpticalDuplicates()));
    }

    /**
     * Method which generates a map of the libraries found tagged in readgroups from the header so they can be serialized as indexes to save space
     */
    public static Map<String, Byte> constructLibraryIndex(final SAMFileHeader header) {
        final List<String> discoveredLibraries = header.getReadGroups().stream()
                .map(r -> { String library = r.getLibrary();
                            return library==null? LibraryIdGenerator.UNKNOWN_LIBRARY : library;} )
                .distinct()
                .collect(Collectors.toList());
        if (discoveredLibraries.size() > 255) {
            throw new GATKException("Detected too many read libraries among read groups header, currently MarkDuplicatesSpark only supports up to 256 unique readgroup libraries but " + discoveredLibraries.size() + " were found");
        }
        final Iterator<Byte> iterator = IntStream.range(0, discoveredLibraries.size()).boxed().map(Integer::byteValue).iterator();
        return Maps.uniqueIndex(iterator, idx -> discoveredLibraries.get(idx));
    }

    /**
     * Method which generates a map of the readgroups from the header so they can be serialized as indexes
     */
    private static Map<String, Short> getHeaderReadGroupIndexMap(final SAMFileHeader header) {
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (readGroups.size() > 65535) {
            throw new GATKException("Detected too many read groups in the header, currently MarkDuplicatesSpark only supports up to 65535 unique readgroup IDs but " + readGroups.size() + " were found");
        }
        if (readGroups.size()==0) {
            throw new UserException.BadInput("Sam file header missing Read Group fields. MarkDuplicatesSpark currently requires reads to be labeled with read group tags, please add read groups tags to your reads");
        }
        final Iterator<Short> iterator = IntStream.range(0, readGroups.size()).boxed().map(Integer::shortValue).iterator();
        return Maps.uniqueIndex(iterator, idx -> readGroups.get(idx).getId() );
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
            throw new GATKException(String.format("MarkDuplicatesSparkUtils.mark() requires input reads to be queryname sorted or querygrouped, yet the header indicated it was in %s order instead", header.getSortOrder()));
        }
        return keyedReads;
    }

    /**
     * Method which takes an RDD of reads that is guaranteed to have every readname group placed together on the same
     * partition and maps those so a JavaPairRDD with the readname as the key.
     */
    private static JavaPairRDD<String, Iterable<IndexPair<GATKRead>>> spanReadsByKey(final JavaRDD<IndexPair<GATKRead>> reads) {
        JavaPairRDD<String, IndexPair<GATKRead>> nameReadPairs = reads.mapToPair(read -> new Tuple2<>(read.getValue().getName(), read));
        return SparkUtils.spanByKey(nameReadPairs).flatMapToPair(namedRead -> {
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
     * Primary landing point for MarkDuplicateSparkRecords:
     *  - Handles separating out hashed keys into into groups by start position/readgroup
     *  - Further separates out MarkDuplicatesSparkRecord by their record objects
     *  - Farms out to methods which handles each group
     *  - Collects the results and returns an iterator
     */
    private static Iterator<NonDuplicatesSparkSqlRecord> markDuplicateRecords(final Iterable<MarkDuplicatesSparkSqlRecord> pairGroups,
                                                                                     final OpticalDuplicateFinder finder, final boolean markOpticalDups) {
        final List<NonDuplicatesSparkSqlRecord> nonDuplicates = Lists.newArrayList();
        final Map<MarkDuplicatesSparkRecord.Type, List<MarkDuplicatesSparkSqlRecord>> stratifiedByType = splitByType(pairGroups);

        // Each key corresponds to either fragments or paired ends, not a mixture of both.
        final List<MarkDuplicatesSparkSqlRecord> emptyFragments = stratifiedByType.get(MarkDuplicatesSparkRecord.Type.EMPTY_FRAGMENT);
        final List<MarkDuplicatesSparkSqlRecord> fragments = stratifiedByType.get(MarkDuplicatesSparkRecord.Type.FRAGMENT);
        final List<MarkDuplicatesSparkSqlRecord> pairs = stratifiedByType.get(MarkDuplicatesSparkRecord.Type.PAIR);
        final List<MarkDuplicatesSparkSqlRecord> passthroughs = stratifiedByType.get(MarkDuplicatesSparkRecord.Type.PASSTHROUGH);

        //empty MarkDuplicatesSparkRecord signify that a pair has a mate somewhere else
        // If there are any non-fragment placeholders at this site, mark everything as duplicates, otherwise compute the best score
        if (Utils.isNonEmpty(fragments) && !Utils.isNonEmpty(emptyFragments)) {
            final NonDuplicatesSparkSqlRecord bestFragment = handleFragments(fragments, finder);
            nonDuplicates.add(bestFragment);
        }

        if (Utils.isNonEmpty(pairs)) {
            nonDuplicates.addAll(handlePairs(pairs, finder, markOpticalDups));
        }

        if (Utils.isNonEmpty(passthroughs)) {
            nonDuplicates.addAll(handlePassthroughs(passthroughs));
        }

        return nonDuplicates.iterator();
    }

    /**
     * split MarkDuplicatesSparkRecord into groups by their type
     */
    private static Map<MarkDuplicatesSparkRecord.Type, List<MarkDuplicatesSparkSqlRecord>> splitByType(Iterable<MarkDuplicatesSparkSqlRecord> duplicateGroup) {
        final EnumMap<MarkDuplicatesSparkRecord.Type, List<MarkDuplicatesSparkSqlRecord>> byType = new EnumMap<>(MarkDuplicatesSparkRecord.Type.class);
        for(MarkDuplicatesSparkSqlRecord pair: duplicateGroup) {
            byType.compute(MarkDuplicatesSparkRecord.Type.values()[pair.getType()], (key, value) -> {
                if (value == null) {
                    final ArrayList<MarkDuplicatesSparkSqlRecord> pairedEnds = new ArrayList<>();
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

    private static List<NonDuplicatesSparkSqlRecord> handlePassthroughs(List<MarkDuplicatesSparkSqlRecord> passthroughs) {
        // Emit the passthrough reads as non-duplicates.
        return passthroughs.stream()
                .map(pair -> new NonDuplicatesSparkSqlRecord(pair.getPartitionIndex(), pair.getName(), MarkDuplicatesSpark.NO_OPTICAL_MARKER))
                .collect(Collectors.toList());
    }

    private static List<NonDuplicatesSparkSqlRecord> handlePairs(final List<MarkDuplicatesSparkSqlRecord> pairs, final OpticalDuplicateFinder finder, final boolean markOpticalDups) {
        // save ourselves the trouble when there are no optical duplicates to worry about
        if (pairs.size() == 1) {
            return Collections.singletonList(new NonDuplicatesSparkSqlRecord(pairs.get(0).getPartitionIndex(), pairs.get(0).getName(), 0));
        }

        List<NonDuplicatesSparkSqlRecord> output = new ArrayList<>();

        final MarkDuplicatesSparkSqlRecord bestPair = pairs.stream()
                .peek(pair -> finder.addLocationInformation(pair.getName(), pair))
                .max(PAIRED_ENDS_SCORE_COMPARATOR)
                .orElseThrow(() -> new GATKException.ShouldNeverReachHereException("There was no best pair because the stream was empty, but it shouldn't have been empty."));

        // Split by orientation and count duplicates in each group separately.
        final Map<Byte, List<MarkDuplicatesSparkSqlRecord>> groupByOrientation = pairs.stream()
                .collect(Collectors.groupingBy(MarkDuplicatesSparkSqlRecord::getOrientationForOpticalDuplicates));
        final int numOpticalDuplicates;
        if (groupByOrientation.containsKey(ReadEnds.FR) && groupByOrientation.containsKey(ReadEnds.RF)) {
            final List<MarkDuplicatesSparkSqlRecord> peFR = new ArrayList<>(groupByOrientation.get(ReadEnds.FR));
            final List<MarkDuplicatesSparkSqlRecord> peRF = new ArrayList<>(groupByOrientation.get(ReadEnds.RF));
            numOpticalDuplicates = countOpticalDuplicates(finder, peFR, bestPair, markOpticalDups? output : null) + countOpticalDuplicates(finder, peRF, bestPair, markOpticalDups? output : null);
        } else {
            numOpticalDuplicates = countOpticalDuplicates(finder, pairs, bestPair, markOpticalDups? output : null);
        }
        output.add(new NonDuplicatesSparkSqlRecord(bestPair.getPartitionIndex(), bestPair.getName(), numOpticalDuplicates));
        return output;
    }

    private static int countOpticalDuplicates(OpticalDuplicateFinder finder, List<MarkDuplicatesSparkSqlRecord> scored, MarkDuplicatesSparkSqlRecord best, List<NonDuplicatesSparkSqlRecord> opticalDuplicateList) {
        final boolean[] opticalDuplicateFlags = finder.findOpticalDuplicates(scored, best);
        int numOpticalDuplicates = 0;
        for (int i = 0; i < opticalDuplicateFlags.length; i++) {
            if (opticalDuplicateFlags[i]) {
                numOpticalDuplicates++;
                if (opticalDuplicateList != null) {
                    opticalDuplicateList.add(new NonDuplicatesSparkSqlRecord(scored.get(i).getPartitionIndex(), scored.get(i).getName(), MarkDuplicatesSpark.OPTICAL_DUPLICATE_MARKER));
                }
            }
        }
        return numOpticalDuplicates;
    }

    /**
     * If there are fragments with no non-fragments overlapping at a site, select the best one according to PAIRED_ENDS_SCORE_COMPARATOR
     */
    private static NonDuplicatesSparkSqlRecord handleFragments(List<MarkDuplicatesSparkSqlRecord> duplicateFragmentGroup, OpticalDuplicateFinder finder) {
        return duplicateFragmentGroup.stream()
                .peek(f -> finder.addLocationInformation(f.getName(), f))
                .max(PAIRED_ENDS_SCORE_COMPARATOR)
                .map(best -> new NonDuplicatesSparkSqlRecord(best.getPartitionIndex(), best.getName(), -1))
                .orElse(null);
    }

    static JavaPairRDD<String, GATKDuplicationMetrics> generateMetrics(final SAMFileHeader header, final JavaRDD<GATKRead> reads) {
        return reads.mapToPair(read -> {
                    final String library = LibraryIdGenerator.getLibraryName(header, read.getReadGroup());
                    GATKDuplicationMetrics metrics = new GATKDuplicationMetrics();
                    metrics.LIBRARY = library;
                    metrics.updateMetrics(read);
                    // NOTE: we use the SAMRecord transientAttribute field here specifically to prevent the already
                    // serialized read from being parsed again here for performance reasons.
                    if (read.getTransientAttribute(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME)!=null) {
                        // NOTE: there is a safety check above in getReadsGroupedByName()
                        metrics.READ_PAIR_OPTICAL_DUPLICATES +=
                                (int)(read.getTransientAttribute(OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME));
                    }
                    return new Tuple2<>(library, metrics);
                })
                .foldByKey(new GATKDuplicationMetrics(), (metricsSum, m) -> {
                    metricsSum.merge(m);
                    if (!metricsSum.LIBRARY.equals(m.LIBRARY)) {
                        throw new GATKException("Two different libraries encountered while summing metrics: " + metricsSum.LIBRARY
                                + " and " + m.LIBRARY);
                    }
                    return metricsSum;
                })
                .mapValues(metrics -> {
                    final GATKDuplicationMetrics copy = metrics.copy();
                    // Divide these by 2 because they are counted for each read
                    // when they should be counted by pair.
                    copy.READ_PAIRS_EXAMINED = metrics.READ_PAIRS_EXAMINED / 2;
                    copy.READ_PAIR_DUPLICATES = metrics.READ_PAIR_DUPLICATES / 2;

                    copy.calculateDerivedFields();
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
    public static void saveMetricsRDD(final MetricsFile<GATKDuplicationMetrics, Double> result, final SAMFileHeader header, final JavaPairRDD<String, GATKDuplicationMetrics> metricsRDD, final String metricsOutputPath) {
        final LibraryIdGenerator libraryIdGenerator = new LibraryIdGenerator(header);

        final Map<String, GATKDuplicationMetrics> nonEmptyMetricsByLibrary = metricsRDD.collectAsMap();           //Unknown Library
        final Map<String, GATKDuplicationMetrics> emptyMapByLibrary = libraryIdGenerator.getMetricsByLibraryMap();//with null

        final List<String> sortedListOfLibraryNames = new ArrayList<>(Sets.union(emptyMapByLibrary.keySet(), nonEmptyMetricsByLibrary.keySet()));
        sortedListOfLibraryNames.sort(Utils.COMPARE_STRINGS_NULLS_FIRST);
        for (final String library : sortedListOfLibraryNames) {
            //if a non-empty exists, take it, otherwise take from the the empties. This is done to include libraries with zero data in them.
            //But not all libraries are listed in the header (esp in testing data) so we union empty and non-empty
            final GATKDuplicationMetrics metricsToAdd = nonEmptyMetricsByLibrary.containsKey(library) ? nonEmptyMetricsByLibrary.get(library) : emptyMapByLibrary.get(library);
            metricsToAdd.calculateDerivedFields();
            result.addMetric(metricsToAdd);
        }

        if (nonEmptyMetricsByLibrary.size() == 1) {
            result.setHistogram(nonEmptyMetricsByLibrary.values().iterator().next().calculateRoiHistogram());
        }

        MetricsUtils.saveMetrics(result, metricsOutputPath);
    }

    /**
     * Comparator for MarkDuplicatesSparkSqlRecord objects by their attributes and strandedness. This comparator is intended to serve as a tiebreaker
     * for the score comparator.
     *
     * It compares two PhysicalLocation  the orientation of the strand, followed by their physical location attributes,
     * and finally as a final tiebreaker the readname lexicographical order.
     *
     * NOTE: Because the original records were grouped by start position, we know that they must be unique once they hit this
     *       comparator and thus we don't need to worry about further tiebreaking for this method.
     */
    public static final class MarkDuplicatesSparkSqlRecordPhysicalLocationComparator implements Comparator<MarkDuplicatesSparkSqlRecord>, Serializable {
        private static final long serialVersionUID = 1L;

        public static final MarkDuplicatesSparkSqlRecordPhysicalLocationComparator INSTANCE = new MarkDuplicatesSparkSqlRecordPhysicalLocationComparator();
        private MarkDuplicatesSparkSqlRecordPhysicalLocationComparator() { }

        @Override
        public int compare( MarkDuplicatesSparkSqlRecord first, MarkDuplicatesSparkSqlRecord second ) {
            int result = 0;

            //This is done to mimic SAMRecordCoordinateComparator's behavior
            if (first.getRead1ReverseStrand() != second.getRead1ReverseStrand()) {
                return first.getRead1ReverseStrand() ? -1: 1;
            }

            if (first.getTile() != second.getTile()) {
                return first.getTile() - second.getTile();
            }
            if (first.getX() != second.getX()) {
                return first.getX() - second.getX();
            }
            if (first.getY() != second.getY()) {
                return first.getY() - second.getY();
            }

            if ( first.getName() != null && second.getName() != null ) {
                result = first.getName().compareTo(second.getName());
            }
            return result;
        }
    }
}