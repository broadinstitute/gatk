package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryInputMetaData;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import scala.Tuple2;

import java.util.List;

/**
 * This deals with the special case where a contig has exactly two alignments
 * and seemingly has the complete alt haplotype assembled.
 * See criteria in {@link AssemblyContigWithFineTunedAlignments#hasIncompletePictureFromTwoAlignments(AlignmentInterval, AlignmentInterval)}.
 *
 * TODO: 1/19/18 see ticket 4189
 *      Exactly how the returned type in {@link SimpleNovelAdjacencyAndChimericAlignmentEvidence} is treated (trusted, or updated, or re-interpreted),
 *      is to be developed.
 */
public final class SimpleNovelAdjacencyInterpreter {

    public static List<VariantContext> makeInterpretation(final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithSimpleChimera,
                                                          final SvDiscoveryInputMetaData svDiscoveryInputMetaData) {

        final JavaPairRDD<SimpleNovelAdjacencyAndChimericAlignmentEvidence, List<SvType>> narlAndAltSeqAndEvidenceAndTypes =
                SimpleNovelAdjacencyInterpreter
                        .inferTypeFromSingleContigSimpleChimera(contigsWithSimpleChimera, svDiscoveryInputMetaData).cache();

        try {
            final List<NovelAdjacencyAndAltHaplotype> narls = narlAndAltSeqAndEvidenceAndTypes.keys()
                    .map(SimpleNovelAdjacencyAndChimericAlignmentEvidence::getNovelAdjacencyReferenceLocations).collect();
            evaluateNarls(svDiscoveryInputMetaData, narls);

            final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast =
                    svDiscoveryInputMetaData.getReferenceData().getReferenceSequenceDictionaryBroadcast();
            final String sampleId = svDiscoveryInputMetaData.getSampleSpecificData().getSampleId();
            final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast = svDiscoveryInputMetaData.getSampleSpecificData().getCnvCallsBroadcast();
            final List<VariantContext> annotatedSimpleVariants =
                    narlAndAltSeqAndEvidenceAndTypes
                            .flatMap(pair ->
                                    pair._1.toVariantContexts(pair._2, sampleId,
                                                referenceSequenceDictionaryBroadcast.getValue(),
                                                cnvCallsBroadcast == null ? null : cnvCallsBroadcast.getValue())
                                            .iterator()
                            )
                            .collect();

            narlAndAltSeqAndEvidenceAndTypes.unpersist();
            return annotatedSimpleVariants;

        } finally {
            narlAndAltSeqAndEvidenceAndTypes.unpersist();
        }
    }

    private static void evaluateNarls(final SvDiscoveryInputMetaData svDiscoveryInputMetaData,
                                      final List<NovelAdjacencyAndAltHaplotype> narls) {
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast =
                svDiscoveryInputMetaData.getReferenceData().getReferenceSequenceDictionaryBroadcast();
        final List<SVInterval> assembledIntervals = svDiscoveryInputMetaData.getSampleSpecificData().getAssembledIntervals();
        final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsArgumentCollection
                discoverStageArgs = svDiscoveryInputMetaData.getDiscoverStageArgs();
        final Logger toolLogger = svDiscoveryInputMetaData.getToolLogger();
        SvDiscoveryUtils.evaluateIntervalsAndNarls(assembledIntervals, narls,
                referenceSequenceDictionaryBroadcast.getValue(), discoverStageArgs, toolLogger);
    }

    /**
     * Filters input assembly contigs that are not strong enough to support an event,
     * then delegates to {@link BreakpointsInference} to infer the reference locations
     * that bound the bi-path bubble in the graph caused by the event,
     * as well as the alternative path encoded in the contig sequence.
     */
    private static JavaPairRDD<SimpleNovelAdjacencyAndChimericAlignmentEvidence, List<SvType>>
    inferTypeFromSingleContigSimpleChimera(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                           final SvDiscoveryInputMetaData svDiscoveryInputMetaData) {

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputMetaData.getReferenceData().getReferenceSequenceDictionaryBroadcast();
        final Broadcast<ReferenceMultiSparkSource> referenceBroadcast = svDiscoveryInputMetaData.getReferenceData().getReferenceBroadcast();

        return
                assemblyContigs
                        .filter(tig -> SimpleChimera
                                .splitPairStrongEnoughEvidenceForCA(tig.getHeadAlignment(), tig.getTailAlignment()))

                        .mapToPair(tig -> getNovelAdjacencyAndEvidence(tig, referenceSequenceDictionaryBroadcast.getValue()))

                        .groupByKey()       // group the same novel adjacency produced by different contigs together

                        .mapToPair(noveltyAndEvidence -> inferType(noveltyAndEvidence, referenceBroadcast));
    }

    private static Tuple2<NovelAdjacencyAndAltHaplotype, SimpleChimera> getNovelAdjacencyAndEvidence
            (final AssemblyContigWithFineTunedAlignments assemblyContigWithFineTunedAlignments,
             final SAMSequenceDictionary refSeqDict) {
        final SimpleChimera simpleChimera = assemblyContigWithFineTunedAlignments.extractSimpleChimera(refSeqDict);
        final byte[] contigSequence = assemblyContigWithFineTunedAlignments.getContigSequence();

        final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype =
                new NovelAdjacencyAndAltHaplotype(simpleChimera, contigSequence, refSeqDict);
        return new Tuple2<>(novelAdjacencyAndAltHaplotype, simpleChimera);
    }

    private static Tuple2<SimpleNovelAdjacencyAndChimericAlignmentEvidence, List<SvType>> inferType
            (final Tuple2<NovelAdjacencyAndAltHaplotype, Iterable<SimpleChimera>> noveltyAndEvidence,
             final Broadcast<ReferenceMultiSparkSource> referenceBroadcast) {

        final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype = noveltyAndEvidence._1;
        final Iterable<SimpleChimera> simpleChimeras = noveltyAndEvidence._2;
        final List<SvType> inferredTypes =
                novelAdjacencyAndAltHaplotype.toSimpleOrBNDTypes(referenceBroadcast.getValue());
        final SimpleNovelAdjacencyAndChimericAlignmentEvidence simpleNovelAdjacencyAndChimericAlignmentEvidence
                = new SimpleNovelAdjacencyAndChimericAlignmentEvidence(novelAdjacencyAndAltHaplotype, simpleChimeras);
        return new Tuple2<>(simpleNovelAdjacencyAndChimericAlignmentEvidence, inferredTypes);
    }
}
