package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.RDDUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;
import scala.Tuple3;
import scala.Tuple4;

import java.util.*;
import java.util.stream.Collectors;

// TODO: 10/30/17 temporary class later to be merged to AssemblyContigAlignmentSignatureClassifier.java after #3752 is in
final class TempMultipleAlignmentsReclassifier {

    // =================================================================================================================

    /**
     * Reclassify assembly contigs based on alignment fine tuning.
     * @return 4 classes:
     *          1) non-informative contigs who after fine tuning has 0 or 1 good alignment left
     *          2) contigs with more than 2 good alignments but doesn't seem to have picture complete as defined by {@link #hasIncomePictureDueToChromosomeHopping(AlignedContig)}
     *          3) contigs with 2 good alignments and bad alignments encoded as strings
     *          4) contigs with more than 2 good alignments and seemingly have picture complete as defined by {@link #hasIncomePictureDueToChromosomeHopping(AlignedContig)}
     */
    static Tuple4<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>, JavaPairRDD<AlignedContig, List<String>>, JavaPairRDD<AlignedContig, List<String>>>
    reClassifyContigsWithMultipleAlignments(final JavaRDD<AlignedContig> localAssemblyContigs,
                                            final int mapQThresholdInclusive, final int uniqReadLenInclusive) {

        final JavaRDD<Tuple3<Boolean, AlignedContig, List<String>>> map = localAssemblyContigs.map(tig -> {
            final Tuple2<List<AlignmentInterval>, List<AlignmentInterval>> goodAndBadAlignments =
                    fineTuneAlignments(tig.alignmentIntervals, mapQThresholdInclusive, uniqReadLenInclusive);
            final List<AlignmentInterval> goodAlignments = goodAndBadAlignments._1;
            final List<String> insertionMappings =
                    goodAndBadAlignments._2.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList());
            return new Tuple3<>(goodAlignments.size() < 2, // after fine tuning, a contig may have no good alignment left, or only 1
                    new AlignedContig(tig.contigName, tig.contigSequence, goodAlignments, tig.hasEquallyGoodAlnConfigurations),
                    insertionMappings);
        });

        // 1st class
        final Tuple2<JavaRDD<Tuple3<Boolean, AlignedContig, List<String>>>, JavaRDD<Tuple3<Boolean, AlignedContig, List<String>>>> split =
                RDDUtils.split(map, TempMultipleAlignmentsReclassifier::isNonInformative, false);
        final JavaRDD<AlignedContig> garbage = split._1.map(Tuple3::_2);

        // 3rd class
        final Tuple2<JavaRDD<Tuple3<Boolean, AlignedContig, List<String>>>, JavaRDD<Tuple3<Boolean, AlignedContig, List<String>>>> split1 =
                RDDUtils.split(split._2, TempMultipleAlignmentsReclassifier::hasOnly2Alignments, false);
        final JavaPairRDD<AlignedContig, List<String>> twoAlignments =
                split1._1.mapToPair(tuple3 -> new Tuple2<>(tuple3._2(), tuple3._3()));

        // 2nd and 4th classes
        final JavaPairRDD<AlignedContig, List<String>> multipleAlignments =
                split1._2.mapToPair(tuple3 -> new Tuple2<>(tuple3._2(), tuple3._3()));
        final JavaPairRDD<AlignedContig, List<String>> multipleAlignmentsIncompletePicture =
                multipleAlignments.filter(pair -> hasIncomePictureDueToChromosomeHopping(pair._1));
        final JavaPairRDD<AlignedContig, List<String>> multipleAlignmentsCompletePicture =
                multipleAlignments.filter(pair -> !hasIncomePictureDueToChromosomeHopping(pair._1));

        return new Tuple4<>(garbage, multipleAlignmentsIncompletePicture.map(Tuple2::_1),
                twoAlignments, multipleAlignmentsCompletePicture);
    }

    private static boolean hasIncomePictureDueToChromosomeHopping(final AlignedContig contigWithMultipleAlignments) {
        Utils.validateArg(contigWithMultipleAlignments.alignmentIntervals.size() > 2,
                "assumption that input contig has more than 2 alignments is violated.\n" +
                        contigWithMultipleAlignments.toString());

        final AlignmentInterval head = contigWithMultipleAlignments.alignmentIntervals.get(0),
                                tail = contigWithMultipleAlignments.alignmentIntervals.get(contigWithMultipleAlignments.alignmentIntervals.size()-1);

        if (!head.referenceSpan.getContig().equals(tail.referenceSpan.getContig()))
            return true;

        if (head.forwardStrand != tail.forwardStrand)
            return true;

        return refSpanAnomaly(head, tail);
    }

    private static boolean refSpanAnomaly(final AlignmentInterval one, final AlignmentInterval two) {
        final SimpleInterval referenceSpanOne = one.referenceSpan,
                             referenceSpanTwo = two.referenceSpan;
        if (referenceSpanOne.contains(referenceSpanTwo) || referenceSpanTwo.contains(referenceSpanOne))
            return true;

        if (one.forwardStrand != two.forwardStrand) {
            return referenceSpanOne.overlaps(referenceSpanTwo);
        } else {
            if (one.forwardStrand) {
                return referenceSpanOne.getStart() > referenceSpanTwo.getStart() &&
                        referenceSpanOne.getStart() <= referenceSpanTwo.getEnd();
            } else {
                return referenceSpanTwo.getStart() > referenceSpanOne.getStart() &&
                        referenceSpanTwo.getStart() <= referenceSpanOne.getEnd();
            }
        }
    }

    private static boolean isNonInformative(final Tuple3<Boolean, AlignedContig, List<String>> tigInfo) {
        return tigInfo._1();
    }

    private static boolean hasOnly2Alignments(final Tuple3<Boolean, AlignedContig, List<String>> tigInfo) {
        return tigInfo._2().alignmentIntervals.size() == 2;
    }

    /**
     * Preprocess provided {@code originalConfiguration} of a particular contig and return a configuration after operation.
     * Note that "original" is meant to be possibly different from the returned configuration,
     * but DOES NOT mean the alignments of the contig as given by the aligner, i.e. the configuration should be
     * one of the best given by {@link FilterLongReadAlignmentsSAMSpark#pickBestConfigurations(AlignedContig, Set, Double)}.
     *
     * Note that the configuration after this preprocessing may not even have the same number of alignment as in input configuration.
     */
    static Tuple2<List<AlignmentInterval>, List<AlignmentInterval>> fineTuneAlignments(final List<AlignmentInterval> originalConfiguration,
                                                                                       final int mapQThresholdInclusive,
                                                                                       final int uniqReadLenInclusive) {
        Utils.validateArg(originalConfiguration.size() > 2,
                "assumption that input configuration to be fine tuned has more than 2 alignments is violated.\n" +
                        originalConfiguration.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList()));

        // two pass, each focusing on removing the alignments of a contig that offers low uniqueness in one sense:

        // first pass is for removing alignments with low REFERENCE UNIQUENESS, using low mapping quality as the criterion
        final List<AlignmentInterval> selectedAlignments = new ArrayList<>(originalConfiguration.size()),
                                      lowUniquenessMappings = new ArrayList<>(originalConfiguration.size());

        for (final AlignmentInterval alignment : originalConfiguration) {
            if (alignment.mapQual >= mapQThresholdInclusive)
                selectedAlignments.add(alignment);
            else
                lowUniquenessMappings.add(alignment);
        }

        // second pass, the slower one, is to remove alignments offering low READ UNIQUENESS,
        // i.e. with only a very short part of the read being uniquely explained by this particular alignment;
        // the steps are:
        //      search bi-directionally until cannot find overlap any more, subtract from it all overlaps.
        //      This gives unique read region it explains.
        //      If this unique read region is "short": shorter than {@code uniqReadLenInclusive}), drop it.

        // each alignment has an entry of a tuple2, one for max overlap front, one for max overlap rear,
        // max overlap front is a tuple2 registering the index and overlap bases count
        final List<Tuple2<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>>> maxOverlapMap =
                getMaxOverlapPairs(selectedAlignments);

        final List<Integer> idxToRemove = new ArrayList<>(selectedAlignments.size());
        for (int i = 0; i < selectedAlignments.size(); ++i) {
            final AlignmentInterval cur = selectedAlignments.get(i);
            final Tuple2<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>> maxOverlapFrontAndRear = maxOverlapMap.get(i);
            final int maxOverlapFront = Math.max(0, maxOverlapFrontAndRear._1._2);
            final int maxOverlapRead = Math.max(0, maxOverlapFrontAndRear._2._2);

            final int uniqReadSpan = cur.endInAssembledContig - cur.startInAssembledContig + 1 - maxOverlapFront - maxOverlapRead;
            if (uniqReadSpan < uniqReadLenInclusive)
                idxToRemove.add(i);
        }

        if ( idxToRemove.isEmpty() )
            return new Tuple2<>(selectedAlignments, lowUniquenessMappings);

        // removing in reverse order so that iterators are not invalidated if we were to remove from start
        final ListIterator<Integer> rit = idxToRemove.listIterator(idxToRemove.size());
        while (rit.hasPrevious()) {
            selectedAlignments.remove( rit.previous().intValue() );
        }

        return new Tuple2<>(selectedAlignments, lowUniquenessMappings);
    }

    /**
     * Each alignment has an entry of a tuple2, one for max overlap front, one for max overlap rear,
     * max overlap front is a tuple2 registering the index and overlap bases count
     */
    private static List<Tuple2<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>>> getMaxOverlapPairs(final List<AlignmentInterval> selectedAlignments) {

        final List<Tuple2<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>>> maxOverlapMap =
                new ArrayList<>(Collections.nCopies(selectedAlignments.size(), new Tuple2<>(new Tuple2<>(-1, -1), new Tuple2<>(-1, -1))));

        for(int i = 0; i < selectedAlignments.size() - 1; ++i) {
            final AlignmentInterval cur = selectedAlignments.get(i);
            int maxOverlap = -1;
            int maxOverlapIdx = -1;
            for (int j = i + 1; j < selectedAlignments.size(); ++j) {
                final int overlap = AlignmentInterval.overlapOnContig(cur, selectedAlignments.get(j));
                if (overlap > 0) {
                    maxOverlap = Math.max(maxOverlap, overlap);
                    maxOverlapIdx = j;
                }
                else // following ones, as guaranteed by the ordering of alignments in the contig, cannot overlap
                    break;
            }
            if (maxOverlap > 0){
                // first set the max_overlap_rear of the current alignment
                Tuple2<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>> oldValue = maxOverlapMap.get(i);
                maxOverlapMap.set(i, new Tuple2<>(oldValue._1, new Tuple2<>(maxOverlapIdx, maxOverlap))); // maxOverlapIdx cannot be -1 here
                // then conditionally set the max_overlap_front of the maxOverlapIdx-th alignment that maximally overlaps with the current alignment
                oldValue = maxOverlapMap.get(maxOverlapIdx);
                if (oldValue._1._2 < maxOverlap)
                    maxOverlapMap.set(maxOverlapIdx, new Tuple2<>(new Tuple2<>(i, maxOverlap), oldValue._2));
            }
        }
        return maxOverlapMap;
    }
}
