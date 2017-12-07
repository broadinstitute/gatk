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
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


public final class InsDelVariantDetector implements VariantDetectorFromLocalAssemblyContigAlignments {

    @Override
    public void inferSvAndWriteVCF(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                   final SvDiscoveryInputData svDiscoveryInputData){

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final String outputPath = svDiscoveryInputData.outputPath;
        final Logger toolLogger = svDiscoveryInputData.toolLogger;

        // TODO: 11/23/17 take insertion mappings from the input and add them to VC
        final JavaPairRDD<byte[], List<ChimericAlignment>> contigSeqAndChimeras =
                assemblyContigs
                        .map( AssemblyContigWithFineTunedAlignments::getSourceContig )
                        .mapToPair(tig -> convertAlignmentIntervalToChimericAlignment(tig,
                                StructuralVariationDiscoveryArgumentCollection.
                                        DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD,
                                StructuralVariationDiscoveryArgumentCollection.
                                        DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH,
                                referenceSequenceDictionaryBroadcast.getValue()));

        final List<VariantContext> annotatedInsDels = produceVariantsFromSimpleChimeras(contigSeqAndChimeras, svDiscoveryInputData);

        SVVCFWriter.writeVCF(annotatedInsDels, outputPath,
                referenceSequenceDictionaryBroadcast.getValue(), toolLogger);
    }

    /**
     * Very similar to {@link ChimericAlignment#parseOneContig(AlignedContig, SAMSequenceDictionary, boolean, int, int, boolean)}, except that
     * badly mapped (MQ < 60) 1st alignment is no longer skipped.
     */
    private static Tuple2<byte[], List<ChimericAlignment>> convertAlignmentIntervalToChimericAlignment (final AlignedContig contig,
                                                                                                        final int mapQualThresholdInclusive,
                                                                                                        final int minAlignmentBlockSize,
                                                                                                        final SAMSequenceDictionary referenceDictionary) {

        final List<AlignmentInterval> alignmentIntervals = contig.alignmentIntervals;
        final Iterator<AlignmentInterval> iterator = alignmentIntervals.iterator();
        AlignmentInterval current = iterator.next();
        final List<ChimericAlignment> results = new ArrayList<>(alignmentIntervals.size() - 1);
        final List<String> insertionMappings = new ArrayList<>();
        while ( iterator.hasNext() ) {
            final AlignmentInterval next = iterator.next();
            if (ChimericAlignment.nextAlignmentMayBeInsertion(current, next, mapQualThresholdInclusive, minAlignmentBlockSize, true)) {
                if (iterator.hasNext()) {
                    insertionMappings.add(next.toPackedString());
                    continue;
                } else {
                    break;
                }
            }

            final ChimericAlignment ca = new ChimericAlignment(current, next, insertionMappings, contig.contigName, referenceDictionary);
            final boolean validStateForThisPath = ca.isNeitherSimpleTranslocationNorIncompletePicture();
            if ( ! validStateForThisPath )
                throw new GATKException.ShouldNeverReachHereException("Mapped assembled contigs are sent down the wrong path: " +
                        "contig suggesting \"translocation\" is sent down the insert/deletion path.\n" + contig.toString());
            results.add(ca);
        }
        return new Tuple2<>(contig.contigSequence,results);
    }

    public static List<VariantContext> produceVariantsFromSimpleChimeras(final JavaPairRDD<byte[], List<ChimericAlignment>> contigSeqAndChimeras,
                                                                         final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final List<SVInterval> assembledIntervals = svDiscoveryInputData.assembledIntervals;
        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast = svDiscoveryInputData.cnvCallsBroadcast;
        final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs = svDiscoveryInputData.discoverStageArgs;
        final String sampleId = svDiscoveryInputData.sampleId;
        final Logger toolLogger = svDiscoveryInputData.toolLogger;

        final JavaPairRDD<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> narlsAndSources =
                contigSeqAndChimeras
                        .flatMapToPair(tigSeqAndChimeras ->
                                discoverNovelAdjacencyFromChimericAlignments(tigSeqAndChimeras,
                                        referenceSequenceDictionaryBroadcast.getValue()))   // a filter-passing contig's alignments may or may not produce novel adjacency, hence flatmap
                        .groupByKey();                                                      // group the same novel adjacency produced by different contigs together

        narlsAndSources.cache();

        SvDiscoveryUtils.evaluateIntervalsAndNarls(assembledIntervals, narlsAndSources, referenceSequenceDictionaryBroadcast,
                discoverStageArgs, toolLogger);

        List<VariantContext> annotatedVariants =
                narlsAndSources
                        .mapToPair(noveltyAndEvidence -> new Tuple2<>(noveltyAndEvidence._1,
                                new Tuple2<>(inferTypeFromNovelAdjacency(noveltyAndEvidence._1), noveltyAndEvidence._2)))       // type inference based on novel adjacency and evidence alignments
                        .map(noveltyTypeAndEvidence ->
                                annotateVariant(                                                                                // annotate the novel adjacency and inferred type
                                        noveltyTypeAndEvidence._1,
                                        noveltyTypeAndEvidence._2._1,
                                        null,
                                        noveltyTypeAndEvidence._2._2,
                                        referenceBroadcast,
                                        referenceSequenceDictionaryBroadcast,
                                        cnvCallsBroadcast,
                                        sampleId))
                        .collect();

        narlsAndSources.unpersist();

        return annotatedVariants;
    }

    //==================================================================================================================

    private static Iterator<Tuple2<NovelAdjacencyReferenceLocations, ChimericAlignment>>
    discoverNovelAdjacencyFromChimericAlignments(final Tuple2<byte[], List<ChimericAlignment>> tigSeqAndChimeras, final SAMSequenceDictionary referenceDictionary) {
        return Utils.stream(tigSeqAndChimeras._2)
                .map(ca -> new Tuple2<>(new NovelAdjacencyReferenceLocations(ca, tigSeqAndChimeras._1, referenceDictionary), ca)).iterator();
    }

    @VisibleForTesting
    public static SimpleSVType inferTypeFromNovelAdjacency(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

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

    static VariantContext annotateVariant(final NovelAdjacencyReferenceLocations novelAdjacency,
                                          final SvType inferredType,
                                          final byte[] altHaplotypeSeq,
                                          final Iterable<ChimericAlignment> chimericAlignments,
                                          final Broadcast<ReferenceMultiSource> broadcastReference,
                                          final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                          final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls,
                                          final String sampleId)
            throws IOException {
        return AnnotatedVariantProducer
                .produceAnnotatedVcFromInferredTypeAndRefLocations(novelAdjacency.leftJustifiedLeftRefLoc,
                        novelAdjacency.leftJustifiedRightRefLoc.getStart(), novelAdjacency.complication,
                        inferredType, altHaplotypeSeq, chimericAlignments,
                        broadcastReference, broadcastSequenceDictionary, broadcastCNVCalls, sampleId);
    }
}
