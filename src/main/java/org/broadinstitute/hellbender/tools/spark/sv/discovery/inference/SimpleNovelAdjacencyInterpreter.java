package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;
import scala.Tuple3;

import java.util.*;
import java.util.stream.Collectors;

/**
 * This deals with the special case where a contig has exactly two alignments
 * and seemingly has the complete alt haplotype assembled.
 * See criteria in {@link AssemblyContigWithFineTunedAlignments#hasIncompletePictureFromMultipleAlignments()}.
 *
 * TODO: 1/19/18 see ticket 4189
 *      Exactly how the returned type in {@link SimpleNovelAdjacency} is treated (trusted, or updated, or re-interpreted),
 *      is to be developed.
 */
public final class SimpleNovelAdjacencyInterpreter {

    static final int MORE_RELAXED_ALIGNMENT_MIN_LENGTH = 30;
    static final int MORE_RELAXED_ALIGNMENT_MIN_MQ = 20;

    public JavaPairRDD<SimpleNovelAdjacency, List<SvType>> inferTypeFromSingleContigSimpleChimera(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                                                                                  final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;

        final JavaRDD<SimpleNovelAdjacency> simpleNovelAdjacencies =
                getSimpleNovelAdjacencyJavaRDD(assemblyContigs, svDiscoveryInputData);

        return simpleNovelAdjacencies
                        .mapToPair(simpleNovelAdjacency ->
                                new Tuple2<>(simpleNovelAdjacency,
                                        inferSimpleOrBNDTypesFromNovelAdjacency(simpleNovelAdjacency,
                                                referenceBroadcast.getValue(), referenceSequenceDictionaryBroadcast.getValue())));
    }

    private JavaRDD<SimpleNovelAdjacency> getSimpleNovelAdjacencyJavaRDD(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                                                         final SvDiscoveryInputData svDiscoveryInputData) {
        final Logger toolLogger = svDiscoveryInputData.toolLogger;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs = svDiscoveryInputData.discoverStageArgs;
        final List<SVInterval> assembledIntervals = svDiscoveryInputData.assembledIntervals;

        final JavaRDD<SimpleNovelAdjacency> simpleNovelAdjacencies =
                assemblyContigs
                        .filter(tig ->
                                ChimericAlignment.splitPairStrongEnoughEvidenceForCA(tig.getSourceContig().alignmentIntervals.get(0),
                                        tig.getSourceContig().alignmentIntervals.get(1),
                                        MORE_RELAXED_ALIGNMENT_MIN_MQ, MORE_RELAXED_ALIGNMENT_MIN_LENGTH))
                        .mapToPair(tig -> {
                            final ChimericAlignment simpleChimera = ChimericAlignment.extractSimpleChimera(tig,
                                    referenceSequenceDictionaryBroadcast.getValue());
                            final byte[] contigSequence = tig.getSourceContig().contigSequence;
                            final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations =
                                    new NovelAdjacencyReferenceLocations(simpleChimera, contigSequence,
                                            referenceSequenceDictionaryBroadcast.getValue());
                            final byte[] altHaplotypeSeq = null; // TODO: 1/21/18 force all subtypes of NovelAdjacencyReferenceLocations to extract alt haplotype
                            final SimpleNovelAdjacency.NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype =
                                    new SimpleNovelAdjacency.NovelAdjacencyAndInferredAltHaptype(novelAdjacencyReferenceLocations, altHaplotypeSeq);
                            return new Tuple2<>(novelAdjacencyAndInferredAltHaptype, simpleChimera);
                        })
                        .groupByKey()       // group the same novel adjacency produced by different contigs together
                        .map(noveltyAndEvidence ->
                                new SimpleNovelAdjacency(noveltyAndEvidence._1, Lists.newArrayList(noveltyAndEvidence._2)));

        SvDiscoveryUtils.evaluateIntervalsAndNarls(assembledIntervals,
                simpleNovelAdjacencies.map(SimpleNovelAdjacency::getNovelAdjacencyReferenceLocations).collect(),
                referenceSequenceDictionaryBroadcast.getValue(), discoverStageArgs, toolLogger);

        return simpleNovelAdjacencies;
    }

    /**
     * Main method: given simple novel adjacency
     *  (that is, affected reference locations, alt haplotype sequence, and chimeric alignment evidence),
     *  infer type.
     *
     * @return the inferred type could be a single entry for simple variants, or a list of two entries with BND mates.
     */
    static List<SvType> inferSimpleOrBNDTypesFromNovelAdjacency(final SimpleNovelAdjacency simpleNovelAdjacency,
                                                                final ReferenceMultiSource reference, final SAMSequenceDictionary referenceDictionary) {

        // based on characteristic of simple chimera, infer type
        final List<ChimericAlignment> alignmentEvidence = simpleNovelAdjacency.getAlignmentEvidence();

        final List<SvType> inferredType;
        final NovelAdjacencyReferenceLocations novelAdjacency = simpleNovelAdjacency.getNovelAdjacencyReferenceLocations();
        if ( alignmentEvidence.stream().allMatch(ChimericAlignment::isLikelySimpleTranslocation) ) { // all indicate simple translocation
            final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForTranslocSuspect =
                    BreakEndVariantType.TransLocBND.getOrderedMates(novelAdjacency,
                            reference, referenceDictionary);
            inferredType = Arrays.asList(orderedMatesForTranslocSuspect._1, orderedMatesForTranslocSuspect._2);
        } else if ( alignmentEvidence.stream().allMatch(ChimericAlignment::isLikelyInvertedDuplication) ) { // all indicate inverted duplication
            inferredType = Collections.singletonList( new SimpleSVType.DuplicationInverted(novelAdjacency) );
        } else if ( alignmentEvidence.stream().map(ca -> ca.strandSwitch).noneMatch(ss -> ss.equals(StrandSwitch.NO_SWITCH)) ) { // all indicate simple (i.e. no duplicate) strand-switch novel adjacency
            final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForInversionSuspect =
                    BreakEndVariantType.InvSuspectBND.getOrderedMates(novelAdjacency, reference);
            inferredType = Arrays.asList(orderedMatesForInversionSuspect._1, orderedMatesForInversionSuspect._2);
        } else if ( alignmentEvidence.stream().allMatch(ChimericAlignment::isNeitherSimpleTranslocationNorIncompletePicture) &&
                alignmentEvidence.stream().map(ca -> ca.strandSwitch).allMatch(ss -> ss.equals(StrandSwitch.NO_SWITCH)) ){ // all point to simple insertion/deletion/small duplication
            inferredType = Collections.singletonList( inferSimpleTypeFromNovelAdjacency(novelAdjacency) );
        } else {
            throw new GATKException
                    .ShouldNeverReachHereException("novel adjacency has its supporting chimeric alignments showing inconsistent behavior\n" +
                    simpleNovelAdjacency.toString());
        }

        return inferredType;
    }

    @VisibleForTesting
    public static SimpleSVType inferSimpleTypeFromNovelAdjacency(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

        final int start = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd();
        final int end = novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
        final StrandSwitch strandSwitch = novelAdjacencyReferenceLocations.strandSwitch;

        final SimpleSVType type;
        if (strandSwitch == StrandSwitch.NO_SWITCH) { // no strand switch happening, so no inversion
            if (start==end) { // something is inserted
                final boolean hasNoDupSeq = !novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation();
                final boolean hasNoInsertedSeq = novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep().isEmpty();
                if (hasNoDupSeq) {
                    if (hasNoInsertedSeq) {
                        throw new GATKException("Something went wrong in type inference, there's suspected insertion happening but no inserted sequence could be inferred "
                                + novelAdjacencyReferenceLocations.toString());
                    } else {
                        type = new SimpleSVType.Insertion(novelAdjacencyReferenceLocations); // simple insertion (no duplication)
                    }
                } else {
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.DuplicationTandem(novelAdjacencyReferenceLocations); // clean expansion of repeat 1 -> 2, or complex expansion
                    } else {
                        type = new SimpleSVType.DuplicationTandem(novelAdjacencyReferenceLocations); // expansion of 1 repeat on ref to 2 repeats on alt with inserted sequence in between the 2 repeats
                    }
                }
            } else {
                final boolean hasNoDupSeq = !novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation();
                final boolean hasNoInsertedSeq = novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep().isEmpty();
                if (hasNoDupSeq) {
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.Deletion(novelAdjacencyReferenceLocations); // clean deletion
                    } else {
                        type = new SimpleSVType.Deletion(novelAdjacencyReferenceLocations); // scarred deletion
                    }
                } else {
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.Deletion(novelAdjacencyReferenceLocations); // clean contraction of repeat 2 -> 1, or complex contraction
                    } else {
                        throw new GATKException("Something went wrong in novel adjacency interpretation: " +
                                " inferring simple SV type from a novel adjacency between two different reference locations, but annotated with both inserted sequence and duplication, which is NOT simple.\n"
                                + novelAdjacencyReferenceLocations.toString());
                    }
                }
            }
        } else {
            type = new SimpleSVType.Inversion(novelAdjacencyReferenceLocations);
        }

        return type;
    }

    // TODO: 1/21/18 hookup at the right place (right now no variants are using this any way because inverted duplication contigs are filtered out)
    static Iterator<Tuple2<Tuple3<NovelAdjacencyReferenceLocations, SimpleSVType.DuplicationInverted, byte[]>, List<ChimericAlignment>>>
    inferInvDupRange(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<Tuple2<ChimericAlignment, byte[]>>> noveltyAndEvidence) {

        final NovelAdjacencyReferenceLocations novelAdjacency = noveltyAndEvidence._1;
        final SimpleSVType.DuplicationInverted duplicationInverted = new SimpleSVType.DuplicationInverted(novelAdjacency);

        // doing this because the same novel adjacency reference locations might be induced by different (probably only slightly) alt haplotypes, so a single group by NARL is not enough
        final Iterable<Tuple2<ChimericAlignment, byte[]>> chimeraAndContigSeq = noveltyAndEvidence._2;
        final Set<Map.Entry<byte[], List<ChimericAlignment>>> alignmentEvidenceGroupedByAltHaplotypeSequence =
                Utils.stream(chimeraAndContigSeq)
                        .collect(
                                Collectors.groupingBy(caAndSeq ->
                                                novelAdjacency.complication.extractAltHaplotypeForInvDup(caAndSeq._1, caAndSeq._2),
                                        Collectors.mapping(caAndSeq -> caAndSeq._1, Collectors.toList())
                                )
                        )
                        .entrySet();

        return alignmentEvidenceGroupedByAltHaplotypeSequence.stream()
                .map(entry -> new Tuple2<>(new Tuple3<>(novelAdjacency, duplicationInverted, entry.getKey()),
                        entry.getValue()))
                .iterator();
    }

}
