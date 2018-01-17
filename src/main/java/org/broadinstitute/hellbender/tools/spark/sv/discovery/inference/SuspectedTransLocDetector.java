package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakEndVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryInputData;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import scala.Tuple2;

public final class SuspectedTransLocDetector implements VariantDetectorFromLocalAssemblyContigAlignments {


    @Override
    public void inferSvAndWriteVCF(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                   final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final String sampleId = svDiscoveryInputData.sampleId;
        final String outputPath = svDiscoveryInputData.outputPath;
        final Logger toolLogger = svDiscoveryInputData.toolLogger;

        toolLogger.info(assemblyContigs.count() + " chimeras indicating strand-switch-less breakpoints");

        // TODO: 11/23/17 take insertion mappings from the input and add them to VC
        final JavaPairRDD<ChimericAlignment, byte[]> chimeraAndSequence =
                assemblyContigs
                        .filter(decoratedTig ->
                                SimpleStrandSwitchVariantDetector.splitPairStrongEnoughEvidenceForCA(
                                                decoratedTig.getSourceContig().alignmentIntervals.get(0), decoratedTig.getSourceContig().alignmentIntervals.get(1),
                                        SimpleStrandSwitchVariantDetector.MORE_RELAXED_ALIGNMENT_MIN_MQ,
                                        0))
                        .mapToPair(decoratedTig ->
                                convertAlignmentIntervalsToChimericAlignment(decoratedTig.getSourceContig(),
                                        referenceSequenceDictionaryBroadcast.getValue())).cache();

        final JavaRDD<VariantContext> annotatedBNDs =
                chimeraAndSequence
                        .mapToPair(pair -> new Tuple2<>(new NovelAdjacencyReferenceLocations(pair._1, pair._2,
                                referenceSequenceDictionaryBroadcast.getValue()), pair._1))
                        .groupByKey()
                        .mapToPair(noveltyAndEvidence -> inferBNDType(noveltyAndEvidence,
                                referenceBroadcast.getValue(),
                                referenceSequenceDictionaryBroadcast.getValue()))
                        .flatMap(noveltyTypeAndEvidence ->
                                AnnotatedVariantProducer
                                        .produceAnnotatedBNDmatesVcFromNovelAdjacency(noveltyTypeAndEvidence._1,
                                                noveltyTypeAndEvidence._2._1, noveltyTypeAndEvidence._2._2,
                                                referenceBroadcast,
                                                referenceSequenceDictionaryBroadcast, sampleId).iterator());

        SVVCFWriter.writeVCF(annotatedBNDs.collect(), outputPath.replace(".vcf", "_transBND.vcf"),
                referenceSequenceDictionaryBroadcast.getValue(), toolLogger);
    }

    private static Tuple2<ChimericAlignment, byte[]> convertAlignmentIntervalsToChimericAlignment(final AlignedContig contig,
                                                                                                  final SAMSequenceDictionary referenceDictionary) {
        return new Tuple2<>(new ChimericAlignment(contig.alignmentIntervals.get(0), contig.alignmentIntervals.get(1),
                EMPTY_INSERTION_MAPPINGS, contig.contigName, referenceDictionary), contig.contigSequence);
    }

    private static Tuple2<NovelAdjacencyReferenceLocations, Tuple2<Tuple2<BreakEndVariantType, BreakEndVariantType>, Iterable<ChimericAlignment>>>
    inferBNDType(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> noveltyAndEvidence,
                 final ReferenceMultiSource reference, final SAMSequenceDictionary referenceDictionary) {

        final NovelAdjacencyReferenceLocations novelAdjacency = noveltyAndEvidence._1;
        final Iterable<ChimericAlignment> chimericAlignments = noveltyAndEvidence._2;
        return new Tuple2<>(novelAdjacency,
                new Tuple2<>(BreakEndVariantType.TransLocBND.getOrderedMates(novelAdjacency, reference, referenceDictionary),
                        chimericAlignments));
    }
}
