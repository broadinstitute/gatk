package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import com.google.common.collect.*;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates.PairedEnds;
import org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates.ReadsKey;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
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
    //Bases below this quality will not be included in picking the best read from a set of duplicates.
    public static final int MIN_BASE_QUAL = 15;

    // Used to set an attribute on the GATKRead marking this read as an optical duplicate.
    private static final String OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME = "OD";

    /**
     * Takes the reads,
     * group them by library, contig, position and orientation,
     * within each group
     *   (a) if there are only fragments, mark all but the highest scoring as duplicates, or,
     *   (b) if at least one is marked as paired, mark all fragments as duplicates.
     *  Note: Emit only the fragments, as the paired reads are handled separately.
     */
    static JavaRDD<GATKRead> transformFragments(final SAMFileHeader header, final JavaRDD<GATKRead> fragments) {
        //
        // Groups reads by keys - keys are tuples of (library, contig, position, orientation).
        //

        JavaPairRDD<String, Iterable<GATKRead>> groupedReads = fragments.mapToPair(read -> {
            read.setIsDuplicate(false);
            return new Tuple2<>(ReadsKey.keyForFragment(header, read), read);
        }).groupByKey();

        return groupedReads.flatMap(v1 -> {
            List<GATKRead> reads = Lists.newArrayList();
            Iterable<GATKRead> readsCopy = Iterables.transform(v1._2(), GATKRead::copy);
            final Map<Boolean, List<GATKRead>> byPairing = StreamSupport.stream(readsCopy.spliterator(), false).collect(Collectors.partitioningBy(
                    read -> ReadUtils.readHasMappedMate(read)
            ));
            // Note the we emit only fragments from this mapper.
            if (byPairing.get(true).isEmpty()) {
                // There are no paired reads, mark all but the highest scoring fragment as duplicate.
                final List<GATKRead> frags = Ordering.natural().reverse().onResultOf((GATKRead read) -> score(read)).immutableSortedCopy(byPairing.get(false));
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
        });
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
    static JavaRDD<GATKRead> transformReads(final SAMFileHeader header, final OpticalDuplicateFinder finder, final JavaRDD<GATKRead> pairs) {

        JavaPairRDD<String, Iterable<GATKRead>> keyedReads =
                pairs.filter(ReadUtils::readHasMappedMate).mapToPair(read -> new Tuple2<>(ReadsKey.keyForRead(header, read), read)).groupByKey();

        JavaPairRDD<String, Iterable<PairedEnds>> keyedPairs = keyedReads.flatMapToPair(stringIterableTuple2 -> {
            List<Tuple2<String, PairedEnds>> out = Lists.newArrayList();
            final List<GATKRead> sorted = Lists.newArrayList(stringIterableTuple2._2());
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

        return markPairedEnds(keyedPairs, finder);
    }

    static JavaRDD<GATKRead> markPairedEnds(final JavaPairRDD<String, Iterable<PairedEnds>> keyedPairs,
                                            final OpticalDuplicateFinder finder) {
        return keyedPairs.flatMap(stringIterableTuple2 -> {
            List<GATKRead> out = Lists.newArrayList();

            Iterable<PairedEnds> pairedEnds = stringIterableTuple2._2();
            final ImmutableListMultimap<Boolean, PairedEnds> paired = Multimaps.index(pairedEnds, pair -> pair.second() != null);

            // As in Picard, unpaired ends left alone.
            for (final PairedEnds pair : paired.get(false)) {
                out.add(pair.first());
            }

            // order by score
            List<PairedEnds> scored = Ordering.natural().reverse().onResultOf((PairedEnds pair) -> pair.score()).sortedCopy(paired.get(true));

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