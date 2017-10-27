package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.utils.RDDUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.EnumMap;
import java.util.List;

final class AssemblyContigAlignmentSignatureClassifier {

    enum RawTypes {
        InsDel,                 // 2 alignments, indicating ins/del, simple duplication expansion/contractions
        IntraChrStrandSwitch,   // breakpoint for events involving intra-chromosome strand switch, particularly for inversion breakpoint suspects
        MappedInsertionBkpt,    // breakpoint for events involving suspected insertions where the inserted sequence could be mapped to either an 1) same-chromosome location without strand switch (tandem duplication suspect) OR 2) diff-chromosome location (MEI/dispersed dup suspect with or without strand switch)
        Cpx,                    // complex sv that seemingly have the complete event assembled
        Incomplete,             // alignment signature indicates the complete picture is not inferable from alignment signature alone
        Ambiguous;              // assembly contig has more than 1 best alignment configuration hence painting multiple pictures
    }

    @VisibleForTesting
    static EnumMap<RawTypes, JavaRDD<AlignedContig>> classifyContigs(final JavaRDD<AlignedContig> contigs,
                                                                     final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                                                     final Logger toolLogger) {

        final EnumMap<RawTypes, JavaRDD<AlignedContig>> contigsByRawTypes = new EnumMap<>(RawTypes.class);

        // long reads with more than 1 best configurations
        contigsByRawTypes.put(RawTypes.Ambiguous,
                contigs.filter(tig -> tig.hasEquallyGoodAlnConfigurations));

        // long reads with only 1 best configuration
        final JavaRDD<AlignedContig> contigsWithOnlyOneBestConfig =
                contigs.filter(lr -> !lr.hasEquallyGoodAlnConfigurations).cache();

        // divert away those likely suggesting cpx sv
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> twoAlignmentAndManyAlignmentTigSets =
                RDDUtils.split(contigsWithOnlyOneBestConfig, AlignedContig::hasOnly2Alignments, true);

        // TODO: 10/29/17 multiple alignments decision tree return contigs with only 2 fine-tuned alignments
        //       the implementation should send back a Tuple2 of {fine-tuned alignments going to be used for analysis, and List(string, ...) of insertion mappings}

        contigsByRawTypes.put(RawTypes.Cpx, twoAlignmentAndManyAlignmentTigSets._2);

        processContigsWithTwoAlignments(twoAlignmentAndManyAlignmentTigSets._1, contigsByRawTypes, broadcastSequenceDictionary);

        return contigsByRawTypes;
    }

    //==================================================================================================================

    private static void processContigsWithTwoAlignments(final JavaRDD<AlignedContig> contigsWithOnlyOneBestConfigAnd2AI,
                                                        final EnumMap<RawTypes, JavaRDD<AlignedContig>> contigsByRawTypes,
                                                        final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary) {
        // first remove overlap
        final JavaRDD<AlignedContig> preprocessedContigs = contigsWithOnlyOneBestConfigAnd2AI.map(
                contig -> {
                    final List<AlignmentInterval> deOverlappedAlignments =
                            ContigAlignmentsModifier.removeOverlap(contig.alignmentIntervals.get(0),
                                    contig.alignmentIntervals.get(1), broadcastSequenceDictionary.getValue());
                    return new AlignedContig(contig.contigName, contig.contigSequence, deOverlappedAlignments,
                            contig.hasEquallyGoodAlnConfigurations);
                }
        );

        // split between the case where both alignments has unique ref span or not
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> hasFullyContainedRefSpanOrNot =
                RDDUtils.split(preprocessedContigs, AssemblyContigAlignmentSignatureClassifier::hasIncompletePicture, false);
        contigsByRawTypes.put(RawTypes.Incomplete, hasFullyContainedRefSpanOrNot._1);

        // split between same chromosome mapping or not
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> sameChrOrNot =
                RDDUtils.split(hasFullyContainedRefSpanOrNot._2, AssemblyContigAlignmentSignatureClassifier::isSameChromosomeMapping, false);
        final JavaRDD<AlignedContig> diffChrBreakpoints = sameChrOrNot._2;

        // split between strand switch or not (NOTE BOTH SAME CHROMOSOME MAPPING)
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> strandSwitchOrNot =
                RDDUtils.split(sameChrOrNot._1, AssemblyContigAlignmentSignatureClassifier::indicatesIntraChrStrandSwitchBkpts, false);
        contigsByRawTypes.put(RawTypes.IntraChrStrandSwitch, strandSwitchOrNot._1);

        // split between ref block switch or not (NOTE BOTH SAME CHROMOSOME MAPPING AND NO STRAND SWITCH)
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> tandemDupBkptOrSimpleInsDel =
                RDDUtils.split(strandSwitchOrNot._2, AssemblyContigAlignmentSignatureClassifier::indicatesIntraChrTandemDupBkpts, false);
        contigsByRawTypes.put(RawTypes.InsDel, tandemDupBkptOrSimpleInsDel._2);

        final JavaRDD<AlignedContig> mappedInsertionBreakpointSuspects = diffChrBreakpoints.union(tandemDupBkptOrSimpleInsDel._1);
        contigsByRawTypes.replace(RawTypes.MappedInsertionBkpt, mappedInsertionBreakpointSuspects);
    }

    //==================================================================================================================

    // TODO: 10/27/17 add tests for these predicates

    static boolean hasIncompletePicture(final AlignedContig contigWithOnlyOneConfigAnd2Aln) {
        return hasIncompletePictureDueToRefSpanContainment(contigWithOnlyOneConfigAnd2Aln)
                ||
                (isSameChromosomeMapping(contigWithOnlyOneConfigAnd2Aln) && hasIncompletePictureDueToOverlappingRefOrderSwitch(contigWithOnlyOneConfigAnd2Aln));
    }

    @VisibleForTesting
    static boolean hasIncompletePictureDueToRefSpanContainment(final AlignedContig contigWithOnlyOneConfigAnd2Aln) {

        final AlignmentInterval one = contigWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(0),
                                two = contigWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(1);
        return one.referenceSpan.contains(two.referenceSpan)
                ||
                two.referenceSpan.contains(one.referenceSpan);
    }

    @VisibleForTesting
    static boolean hasIncompletePictureDueToOverlappingRefOrderSwitch(final AlignedContig contigWithOnlyOneConfigAnd2AlnToSameChr) {

        final AlignmentInterval one = contigWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(0),
                                two = contigWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(1);
        final SimpleInterval referenceSpanOne = one.referenceSpan,
                             referenceSpanTwo = two.referenceSpan;

        if (referenceSpanOne.contains(referenceSpanTwo) || referenceSpanTwo.contains(referenceSpanOne))
            return true;

        // TODO: 10/29/17 this obsoletes the inverted duplication call code we have now, but those could be used to figure out how to annotate which known ref regions are invert duplicated 
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

    /**
     * @param contigWithOnlyOneConfigAnd2Aln assumed to have only two alignments
     */
    @VisibleForTesting
    static boolean isSameChromosomeMapping(final AlignedContig contigWithOnlyOneConfigAnd2Aln) {
        Utils.validateArg(!hasIncompletePictureDueToRefSpanContainment(contigWithOnlyOneConfigAnd2Aln),
                "assumption that input contig has alignments whose ref span is not completely covered by the other's is violated \n" +
                        contigWithOnlyOneConfigAnd2Aln.toString());
        return contigWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(0).referenceSpan.getContig()
                .equals(contigWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(1).referenceSpan.getContig());
    }

    /**
     * @param contigWithOnlyOneConfigAnd2AlnToSameChr assumed to have only two alignments mapped to the same chromosome
     */
    @VisibleForTesting
    static boolean indicatesIntraChrStrandSwitchBkpts(final AlignedContig contigWithOnlyOneConfigAnd2AlnToSameChr) {
        Utils.validateArg(isSameChromosomeMapping(contigWithOnlyOneConfigAnd2AlnToSameChr),
                "assumption that input contig's 2 alignments map to the same chr is violated. \n" +
                        contigWithOnlyOneConfigAnd2AlnToSameChr.toString());
        return contigWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(0).forwardStrand
                ^
                contigWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(1).forwardStrand;
    }

    /**
     * @param contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch assumed to have only two alignments mapped to the same chromosome and the same strand
     */
    @VisibleForTesting
    static boolean indicatesIntraChrTandemDupBkpts(final AlignedContig contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch) {

        final AlignmentInterval one = contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch.alignmentIntervals.get(0),
                                two = contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch.alignmentIntervals.get(1);
        return ChimericAlignment.isLikelySimpleTranslocation(one, two, StrandSwitch.NO_SWITCH);
    }

    //==================================================================================================================

    // TODO: 11/17/17 salvation on assembly contigs that 1) has ambiguous "best" configuration, and 2) has incomplete picture; and flag accordingly

}
