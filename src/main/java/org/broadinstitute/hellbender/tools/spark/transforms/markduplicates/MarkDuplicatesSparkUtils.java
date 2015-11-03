package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimaps;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.DuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.LibraryIdGenerator;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
import org.broadinstitute.hellbender.utils.read.markduplicates.PairedEnds;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;
import scala.Tuple2;

import java.io.Serializable;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Utility classes and functions for Mark Duplicates.
 */
public class MarkDuplicatesSparkUtils {

    // Used to set an attribute on the GATKRead marking this read as an optical duplicate.
    private static final String OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME = "OD";

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
    static JavaRDD<GATKRead> transformReads(final SAMFileHeader header, final OpticalDuplicateFinder finder, final JavaRDD<GATKRead> reads, final int parallelism) {

        JavaPairRDD<String, Iterable<GATKRead>> keyedReads =
                reads.mapToPair(read -> new Tuple2<>(ReadsKey.keyForRead(header, read), read)).groupByKey(parallelism);

        JavaPairRDD<String, Iterable<PairedEnds>> keyedPairs = keyedReads.flatMapToPair(keyedRead -> {
            List<Tuple2<String, PairedEnds>> out = Lists.newArrayList();
            // Write each read out as a pair with only the first slot filled
            for (GATKRead read : keyedRead._2()) {
                read.setIsDuplicate(false);
                final PairedEnds pair = PairedEnds.of(read);
                out.add(new Tuple2<>(pair.keyForFragment(header), pair));
            }
            // Write each paired read with a mapped mate as a pair
            final List<GATKRead> sorted = Lists.newArrayList(Iterables.filter(keyedRead._2(), read -> ReadUtils.readHasMappedMate(read)));
            sorted.sort(new GATKOrder(header));
            PairedEnds pair = null;
            //Records are sorted, we iterate over them and pair them up.
            for (final GATKRead record : sorted) {
                if (pair == null) {                                //first in pair
                    pair = PairedEnds.of(record);
                } else {                                           //second in pair
                    pair.and(record);
                    out.add(new Tuple2<>(pair.key(header), pair));
                    pair = null;                                   //back to first
                }
            }
            if (pair != null) {                                    //left over read
                out.add(new Tuple2<>(pair.key(header), pair));
            }
            return out;
        }).groupByKey();

        return markPairedEnds(keyedPairs, finder, header);
    }

    static JavaRDD<GATKRead> markPairedEnds(final JavaPairRDD<String, Iterable<PairedEnds>> keyedPairs,
                                            final OpticalDuplicateFinder finder, final SAMFileHeader header) {
        return keyedPairs.flatMap(keyedPair -> {
            Iterable<PairedEnds> pairedEnds = keyedPair._2();
            final ImmutableListMultimap<Boolean, PairedEnds> paired = Multimaps.index(pairedEnds, pair -> pair.second() != null);

            // Each key corresponds to either fragments or paired ends, not a mixture of both.

            if (ReadsKey.isFragment(keyedPair._1())) { // fragments
                return handleFragments(pairedEnds, header);
            }

            List<GATKRead> out = Lists.newArrayList();

            // As in Picard, unpaired ends left alone.
            for (final PairedEnds pair : paired.get(false)) {
                out.add(pair.first());
            }

            // Order by score using ReadCoordinateComparator for tie-breaking.
            Comparator<PairedEnds> pairedEndsComparator =
                    Comparator.comparing(PairedEnds::score).reversed()
                            .thenComparing((o1, o2) -> new ReadCoordinateComparator(header).compare(o1.first(), o2.first()));            List <PairedEnds> scored = paired.get(true).stream().sorted(pairedEndsComparator).collect(Collectors.toList());

            final PairedEnds best = Iterables.getFirst(scored, null);
            if (best == null) {
                return out;
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
                out.add(pair.first());
                out.add(pair.second());
            }
            return out;
        });
    }

    private static List<GATKRead> handleFragments(Iterable<PairedEnds> pairedEnds, final SAMFileHeader header) {
        List<GATKRead> reads = Lists.newArrayList();

        final Iterable<GATKRead> transform = Iterables.transform(pairedEnds, pair -> pair.first());
        Iterable<GATKRead> readsCopy = Iterables.transform(transform, GATKRead::copy);
        final Map<Boolean, List<GATKRead>> byPairing = StreamSupport.stream(readsCopy.spliterator(), false).collect(Collectors.partitioningBy(
                read -> ReadUtils.readHasMappedMate(read)
        ));
        // Note the we emit only fragments from this mapper.
        if (byPairing.get(true).isEmpty()) {
            // There are no paired reads, mark all but the highest scoring fragment as duplicate.
            Comparator<GATKRead> fragmentsComparator = Comparator.comparing((java.util.function.Function<GATKRead, Integer>) (r) -> -MarkDuplicatesUtils.scoreForRead(r)).thenComparing(new ReadCoordinateComparator(header));
            List <GATKRead> frags = byPairing.get(false).stream().sorted(fragmentsComparator).collect(Collectors.toList());
            if (!frags.isEmpty()) {
                reads.add(frags.get(0));                        //highest score - just emit
                for (final GATKRead record : Iterables.skip(frags, 1)) {  //lower   scores - mark as dups and emit
                    record.setIsDuplicate(true);
                    reads.add(record);
                }
            }
        } else {
            // There are paired ends so we mark all fragments as duplicates.
            for (final GATKRead record : byPairing.get(false)) {
                record.setIsDuplicate(true);
                reads.add(record);
            }
        }
        return reads;
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

    public static void saveMetricsRDD(final JavaPairRDD<String, DuplicationMetrics> metricsRDD, final String metricsOutputPath, AuthHolder authHolder) {
        final MetricsFile<DuplicationMetrics, Double> metrics = new MetricsFile<>();
        for (final Map.Entry<String, DuplicationMetrics> entry : metricsRDD.collectAsMap().entrySet()) {
            metrics.addMetric(entry.getValue());
        }

        MetricsUtils.saveMetrics(metrics, metricsOutputPath, authHolder );
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