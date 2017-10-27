package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.utils.RDDUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;
import scala.Tuple4;

import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumMap;
import java.util.List;

final class AssemblyContigAlignmentSignatureClassifier {

    enum RawTypes {
        InsDel,                 // 2 alignments, indicating ins/del, simple duplication expansion/contractions
        IntraChrStrandSwitch,   // intra-chromosome strand switch breakpoints
        TandemDupOrMEIBkpt,     // intra-chromosome tandem duplications without strand switch OR inter-chromosome MEI (with or without strand switch) breakpoints
        Cpx,                    // complex sv that seemingly have the complete event assembled
        Incomplete,             // alignment signature indicates the complete picture is not inferable from alignment signature alone
        Ambiguous,              // assembly contig has more than 1 best alignment configuration hence painting multiple pictures
        MisAssemblySuspect;     // suspected to be misassembly due to alignment signature that despite multiple alignments, no or only 1 good alignment
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
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> hasTwoAlignmentsOrMore =
                RDDUtils.split(contigsWithOnlyOneBestConfig, AssemblyContigAlignmentSignatureClassifier::hasOnly2Alignments, true);

        final Tuple4<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>, JavaPairRDD<AlignedContig, List<String>>, JavaPairRDD<AlignedContig, List<String>>> fineTunedContigs =
                TempMultipleAlignmentsReclassifier.reClassifyContigsWithMultipleAlignments(hasTwoAlignmentsOrMore._2,
                        +                20, 10);
        contigsByRawTypes.put(RawTypes.MisAssemblySuspect, fineTunedContigs._1());
        contigsByRawTypes.put(RawTypes.Incomplete, fineTunedContigs._2());
        contigsByRawTypes.put(RawTypes.Cpx, fineTunedContigs._4().keys());

        final List<String> emptyInsertionMappingStrings = new ArrayList<>();
        final JavaPairRDD<AlignedContig, List<String>> contigsWith2AlignmentsWithOptionallyKnownInsertionMappings =
                hasTwoAlignmentsOrMore._1.mapToPair(tig -> new Tuple2<>(tig, emptyInsertionMappingStrings))
                        .union(fineTunedContigs._3());
        decisionTreeForContigsWithTwoAlignments(contigsWith2AlignmentsWithOptionallyKnownInsertionMappings, contigsByRawTypes, broadcastSequenceDictionary);

        return contigsByRawTypes;
    }

    //==================================================================================================================

    // needs this due to Java naive lambda NOT Spark serializable, and we need a serializable callable
    private static boolean hasOnly2Alignments(final AlignedContig contigWithOnlyOneConfig) {
        return contigWithOnlyOneConfig.alignmentIntervals.size() == 2;
    }

    private static void decisionTreeForContigsWithTwoAlignments(final JavaPairRDD<AlignedContig, List<String>> contigsWithOnlyOneBestConfigAnd2AIWithOptionallyKnownInsertionMappings,
                                                                final EnumMap<RawTypes, JavaRDD<AlignedContig>> contigsByRawTypes,
                                                                final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary) {
        // first remove overlap
        final JavaRDD<AlignedContig> preprocessedContigs =
                contigsWithOnlyOneBestConfigAnd2AIWithOptionallyKnownInsertionMappings.keys().map(
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
        final JavaRDD<AlignedContig> incompletePicture = contigsByRawTypes.get(RawTypes.Incomplete).union(hasFullyContainedRefSpanOrNot._1);
        contigsByRawTypes.put(RawTypes.Incomplete, incompletePicture);

        // split between same chromosome mapping or not
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> sameChrOrNot =
                RDDUtils.split(hasFullyContainedRefSpanOrNot._2, AssemblyContigAlignmentSignatureClassifier::isSameChromosomeMapping, false);
        contigsByRawTypes.put(RawTypes.TandemDupOrMEIBkpt, sameChrOrNot._2);

        // split between strand switch or not (NOTE BOTH SAME CHROMOSOME MAPPING)
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> strandSwitchOrNot =
                RDDUtils.split(sameChrOrNot._1, AssemblyContigAlignmentSignatureClassifier::indicatesIntraChrStrandSwitchBkpts, false);
        contigsByRawTypes.put(RawTypes.IntraChrStrandSwitch, strandSwitchOrNot._1);

        // split between ref block switch or not (NOTE BOTH SAME CHROMOSOME MAPPING AND NO STRAND SWITCH)
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> tandemDupBkptOrSimpleInsDel =
                RDDUtils.split(strandSwitchOrNot._2, AssemblyContigAlignmentSignatureClassifier::indicatesIntraChrTandemDupBkpts, false);
        contigsByRawTypes.put(RawTypes.InsDel, tandemDupBkptOrSimpleInsDel._2);
        contigsByRawTypes.replace(RawTypes.TandemDupOrMEIBkpt, contigsByRawTypes.get(RawTypes.TandemDupOrMEIBkpt).union(tandemDupBkptOrSimpleInsDel._1));
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
        Utils.validateArg(hasOnly2Alignments(contigWithOnlyOneConfigAnd2Aln),
                "assumption that input contig has only 2 alignments is violated. \n" +
                        contigWithOnlyOneConfigAnd2Aln.toString());
        final AlignmentInterval one = contigWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(0),
                                two = contigWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(1);
        return one.referenceSpan.contains(two.referenceSpan)
                ||
                two.referenceSpan.contains(one.referenceSpan);
    }

    @VisibleForTesting
    static boolean hasIncompletePictureDueToOverlappingRefOrderSwitch(final AlignedContig contigWithOnlyOneConfigAnd2AlnToSameChr) {
        Utils.validateArg(isSameChromosomeMapping(contigWithOnlyOneConfigAnd2AlnToSameChr),
                "assumption that input contig's 2 alignments map to the same chr is violated. \n" +
                        contigWithOnlyOneConfigAnd2AlnToSameChr.toString());

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
        Utils.validateArg( !indicatesIntraChrStrandSwitchBkpts(contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch),
                "assumption that input contig's 2 alignments map to the same chr is violated. \n" +
                        contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch.toString());

        final AlignmentInterval one = contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch.alignmentIntervals.get(0),
                                two = contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch.alignmentIntervals.get(1);
        final ChimericAlignment tempCA = new ChimericAlignment(one, two, dummyInsertionMappingStringList,
                "dummy", null); // passing null because alignments are guanranteed to be same-chr mapping
        return !tempCA.isNotSimpleTranslocation();
    }

    private static final List<String> dummyInsertionMappingStringList = Collections.emptyList();

    //==================================================================================================================

    // TODO: 10/27/17 salvation on these reads, and flag accordingly
    static EnumMap<RawTypes, JavaRDD<AlignedContig>> dumpsterDivingInAmbiguous() {
        return new EnumMap<>(RawTypes.class);
    }

    static void dumpsterDivingInIncompletePicture(final JavaRDD<AlignedContig> incompletePictures) {

    }
}
