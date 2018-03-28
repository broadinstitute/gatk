package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import scala.Tuple2;

import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
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

    static final int MORE_RELAXED_ALIGNMENT_MIN_LENGTH = 30;
    static final int MORE_RELAXED_ALIGNMENT_MIN_MQ = 20;

    public static List<VariantContext> makeInterpretation(final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithSimpleChimera,
                                                          final SvDiscoveryInputMetaData svDiscoveryInputMetaData) {

        final JavaPairRDD<SimpleNovelAdjacencyAndChimericAlignmentEvidence, List<SvType>> narlAndAltSeqAndEvidenceAndTypes =
                SimpleNovelAdjacencyInterpreter
                        .inferTypeFromSingleContigSimpleChimera(contigsWithSimpleChimera, svDiscoveryInputMetaData).cache();

        try {
            final List<NovelAdjacencyAndAltHaplotype> narls = narlAndAltSeqAndEvidenceAndTypes.keys()
                    .map(SimpleNovelAdjacencyAndChimericAlignmentEvidence::getNovelAdjacencyReferenceLocations).collect();
            evaluateNarls(svDiscoveryInputMetaData, narls);

            final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputMetaData.referenceData.referenceBroadcast;
            final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast =
                    svDiscoveryInputMetaData.referenceData.referenceSequenceDictionaryBroadcast;
            final String sampleId = svDiscoveryInputMetaData.sampleSpecificData.sampleId;
            final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast = svDiscoveryInputMetaData.sampleSpecificData.cnvCallsBroadcast;
            final List<VariantContext> annotatedSimpleVariants =
                    narlAndAltSeqAndEvidenceAndTypes
                            .flatMap(pair ->
                                    turnIntoVariantContexts(pair, sampleId, referenceBroadcast,
                                            referenceSequenceDictionaryBroadcast, cnvCallsBroadcast)
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
                svDiscoveryInputMetaData.referenceData.referenceSequenceDictionaryBroadcast;
        final List<SVInterval> assembledIntervals = svDiscoveryInputMetaData.sampleSpecificData.assembledIntervals;
        final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection
                discoverStageArgs = svDiscoveryInputMetaData.discoverStageArgs;
        final Logger toolLogger = svDiscoveryInputMetaData.toolLogger;
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

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputMetaData.referenceData.referenceSequenceDictionaryBroadcast;
        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputMetaData.referenceData.referenceBroadcast;

        return
                assemblyContigs
                        .filter(tig -> SimpleChimera
                                .splitPairStrongEnoughEvidenceForCA(tig.getHeadAlignment(), tig.getTailAlignment(),
                                        MORE_RELAXED_ALIGNMENT_MIN_MQ, MORE_RELAXED_ALIGNMENT_MIN_LENGTH))

                        .mapToPair(tig -> getNovelAdjacencyAndEvidence(tig, referenceSequenceDictionaryBroadcast.getValue()))

                        .groupByKey()       // group the same novel adjacency produced by different contigs together

                        .mapToPair(noveltyAndEvidence -> inferType(noveltyAndEvidence, referenceSequenceDictionaryBroadcast, referenceBroadcast));
    }

    private static Tuple2<NovelAdjacencyAndAltHaplotype, SimpleChimera> getNovelAdjacencyAndEvidence
            (final AssemblyContigWithFineTunedAlignments assemblyContigWithFineTunedAlignments,
             final SAMSequenceDictionary refSeqDict) {
        final SimpleChimera simpleChimera = extractSimpleChimera(assemblyContigWithFineTunedAlignments, refSeqDict);
        final byte[] contigSequence = assemblyContigWithFineTunedAlignments.getContigSequence();

        final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype =
                new NovelAdjacencyAndAltHaplotype(simpleChimera, contigSequence, refSeqDict);
        return new Tuple2<>(novelAdjacencyAndAltHaplotype, simpleChimera);
    }

    private static Tuple2<SimpleNovelAdjacencyAndChimericAlignmentEvidence, List<SvType>> inferType
            (final Tuple2<NovelAdjacencyAndAltHaplotype, Iterable<SimpleChimera>> noveltyAndEvidence,
             final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast,
             final Broadcast<ReferenceMultiSource> referenceBroadcast) {

        final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype = noveltyAndEvidence._1;
        final Iterable<SimpleChimera> simpleChimeras = noveltyAndEvidence._2;
        final List<SvType> inferredTypes =
                novelAdjacencyAndAltHaplotype.toSimpleOrBNDTypes(referenceBroadcast.getValue(),
                        referenceSequenceDictionaryBroadcast.getValue());
        final SimpleNovelAdjacencyAndChimericAlignmentEvidence simpleNovelAdjacencyAndChimericAlignmentEvidence
                = new SimpleNovelAdjacencyAndChimericAlignmentEvidence(novelAdjacencyAndAltHaplotype, simpleChimeras);
        return new Tuple2<>(simpleNovelAdjacencyAndChimericAlignmentEvidence, inferredTypes);
    }

    /**
     * @return a simple chimera indicated by the alignments of the input contig;
     *         if the input chimeric alignments are not strong enough to support an CA, a {@code null} is returned
     *
     * @throws IllegalArgumentException if the input contig doesn't have exactly two good input alignments
     */
    @VisibleForTesting
    static SimpleChimera extractSimpleChimera(final AssemblyContigWithFineTunedAlignments contig,
                                              final SAMSequenceDictionary referenceDictionary) {
        if ( ! contig.hasOnly2GoodAlignments() )
            throw new IllegalArgumentException("assembly contig sent to the wrong path: assumption that contig has only 2 good alignments is violated for\n" +
                    contig.toString());

        final AlignmentInterval alignmentOne = contig.getAlignments().get(0);
        final AlignmentInterval alignmentTwo = contig.getAlignments().get(1);

        return new SimpleChimera(alignmentOne, alignmentTwo, contig.getInsertionMappings(),
                contig.getContigName(), contig.getSAtagForGoodMappingToNonCanonicalChromosome(), referenceDictionary);
    }

    /**
     * This implementation is the 1st step going towards allowing re-interpretation,
     * below we simply take the inferred type and turn it to a VC,
     * future implementation may integrate other types of evidence and re-interpret if necessary
     */
    private static Iterator<VariantContext> turnIntoVariantContexts(final Tuple2<SimpleNovelAdjacencyAndChimericAlignmentEvidence, List<SvType>> pair,
                                                                    final String sampleId,
                                                                    final Broadcast<ReferenceMultiSource> referenceBroadcast,
                                                                    final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast,
                                                                    final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast)
            throws IOException {
        final SimpleNovelAdjacencyAndChimericAlignmentEvidence simpleNovelAdjacencyAndChimericAlignmentEvidence = pair._1;
        final List<SvType> svTypes = pair._2;
        if( svTypes.isEmpty() || svTypes.size() > 2 ) {
            throw new GATKException("Wrong number of variants sent for analysis: " + pair._2.toString() +
                    "\nWe currently only support 1 (symbolic simple or CPX) or 2 (BND mate pairs) variants for producing annotated variants.");
        }
        if ( ! svTypes.get(0).isBreakEndOnly() ) { // simple SV type
            final NovelAdjacencyAndAltHaplotype narl = simpleNovelAdjacencyAndChimericAlignmentEvidence.getNovelAdjacencyReferenceLocations();
            final List<SimpleChimera> contigEvidence = simpleNovelAdjacencyAndChimericAlignmentEvidence.getAlignmentEvidence();

            if ( svTypes.size() == 2 ) { // RPL case with both path >= 50 bp
                final SvType firstVar = svTypes.get(0);
                final SvType secondVar = svTypes.get(1);
                final Tuple2<SvType, SvType> linkedVariants = new Tuple2<>(firstVar, secondVar);
                return AnnotatedVariantProducer.produceAnnotatedAndLinkedVcFromNovelAdjacency(linkedVariants,
                        simpleNovelAdjacencyAndChimericAlignmentEvidence,
                        referenceBroadcast, referenceSequenceDictionaryBroadcast, cnvCallsBroadcast, sampleId,
                        GATKSVVCFConstants.LINK).iterator();
            } else {
                final SvType inferredType = svTypes.get(0);

                final VariantContext variantContext = AnnotatedVariantProducer
                        .produceAnnotatedVcFromInferredTypeAndRefLocations(
                                narl, inferredType, contigEvidence,
                                referenceBroadcast, referenceSequenceDictionaryBroadcast, cnvCallsBroadcast, sampleId);
                return Collections.singletonList(variantContext).iterator();
            }
        } else { // BND mate pair
            final BreakEndVariantType firstMate = (BreakEndVariantType) svTypes.get(0);
            final BreakEndVariantType secondMate = (BreakEndVariantType) svTypes.get(1);

            final Tuple2<SvType, SvType> bndMates = new Tuple2<>(firstMate, secondMate);
            final List<VariantContext> variantContexts = AnnotatedVariantProducer
                    .produceAnnotatedAndLinkedVcFromNovelAdjacency(
                            bndMates, simpleNovelAdjacencyAndChimericAlignmentEvidence,
                            referenceBroadcast, referenceSequenceDictionaryBroadcast, cnvCallsBroadcast, sampleId,
                            GATKSVVCFConstants.BND_MATEID_STR);
            return variantContexts.iterator();
        }
    }
}
