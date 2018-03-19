package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import scala.Tuple2;

import java.util.List;

/**
 * This deals with the special case where a contig has exactly two alignments
 * and seemingly has the complete alt haplotype assembled.
 * See criteria in {@link AssemblyContigWithFineTunedAlignments#hasIncompletePictureFromTwoAlignments()} ()}.
 *
 * TODO: 1/19/18 see ticket 4189
 *      Exactly how the returned type in {@link SimpleNovelAdjacencyAndChimericAlignmentEvidence} is treated (trusted, or updated, or re-interpreted),
 *      is to be developed.
 */
public final class SimpleNovelAdjacencyInterpreter {

    static final int MORE_RELAXED_ALIGNMENT_MIN_LENGTH = 30;
    static final int MORE_RELAXED_ALIGNMENT_MIN_MQ = 20;

    public JavaPairRDD<SimpleNovelAdjacencyAndChimericAlignmentEvidence, List<SvType>>
    inferTypeFromSingleContigSimpleChimera(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                           final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;

        final JavaRDD<SimpleNovelAdjacencyAndChimericAlignmentEvidence> simpleNovelAdjacencies =
                getSimpleNovelAdjacencyAndChimeraEvidence(assemblyContigs, svDiscoveryInputData);

        return simpleNovelAdjacencies
                        .mapToPair(simpleNovelAdjacencyAndChimericAlignmentEvidence ->
                                new Tuple2<>(simpleNovelAdjacencyAndChimericAlignmentEvidence,
                                        simpleNovelAdjacencyAndChimericAlignmentEvidence.getNovelAdjacencyReferenceLocations()
                                                .toSimpleOrBNDTypes(referenceBroadcast.getValue(), referenceSequenceDictionaryBroadcast.getValue())));
    }

    /**
     * Filters input assembly contigs that are not strong enough to support an event,
     * then delegates to {@link BreakpointsInference} to infer the reference locations
     * that bound the bi-path bubble in the graph caused by the event,
     * as well as the alternative path encoded in the contig sequence.
     */
    private JavaRDD<SimpleNovelAdjacencyAndChimericAlignmentEvidence>
    getSimpleNovelAdjacencyAndChimeraEvidence(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                              final SvDiscoveryInputData svDiscoveryInputData) {
        final Logger toolLogger = svDiscoveryInputData.toolLogger;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs = svDiscoveryInputData.discoverStageArgs;
        final List<SVInterval> assembledIntervals = svDiscoveryInputData.assembledIntervals;

        final JavaRDD<SimpleNovelAdjacencyAndChimericAlignmentEvidence> simpleNovelAdjacencies =
                assemblyContigs
                        .filter(tig -> ChimericAlignment
                                .splitPairStrongEnoughEvidenceForCA(
                                        tig.getHeadAlignment(),
                                        tig.getTailAlignment(),
                                        MORE_RELAXED_ALIGNMENT_MIN_MQ, MORE_RELAXED_ALIGNMENT_MIN_LENGTH))
                        .mapToPair(tig -> {
                            final SAMSequenceDictionary refSeqDict = referenceSequenceDictionaryBroadcast.getValue();
                            final ChimericAlignment simpleChimera = ChimericAlignment.extractSimpleChimera(tig, refSeqDict);
                            final byte[] contigSequence = tig.getContigSequence();

                            final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype =
                                    new NovelAdjacencyAndAltHaplotype(simpleChimera, contigSequence, refSeqDict);
                            return new Tuple2<>(novelAdjacencyAndAltHaplotype,
                                    new Tuple2<>(simpleChimera, tig.getSAtagForGoodMappingToNonCanonicalChromosome()));
                        })
                        .groupByKey()       // group the same novel adjacency produced by different contigs together
                        .map(noveltyAndEvidence ->
                                new SimpleNovelAdjacencyAndChimericAlignmentEvidence(noveltyAndEvidence._1, noveltyAndEvidence._2));

        SvDiscoveryUtils.evaluateIntervalsAndNarls(assembledIntervals,
                simpleNovelAdjacencies.map(SimpleNovelAdjacencyAndChimericAlignmentEvidence::getNovelAdjacencyReferenceLocations).collect(),
                referenceSequenceDictionaryBroadcast.getValue(), discoverStageArgs, toolLogger);

        return simpleNovelAdjacencies;
    }
}
