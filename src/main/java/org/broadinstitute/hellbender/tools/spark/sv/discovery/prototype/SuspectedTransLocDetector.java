package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import scala.Tuple2;

import java.util.Collections;
import java.util.List;

final class SuspectedTransLocDetector implements VariantDetectorFromLocalAssemblyContigAlignments {

    @SuppressWarnings("unchecked")
    private static final List<String> EMPTY_INSERTION_MAPPINGS = Collections.EMPTY_LIST;

    @Override
    public void inferSvAndWriteVCF(final String vcfOutputFileName, final String sampleId, final JavaRDD<AlignedContig> localAssemblyContigs,
                                   final Broadcast<ReferenceMultiSource> broadcastReference,
                                   final Broadcast<SAMSequenceDictionary> refSequenceDictionary,
                                   final Logger toolLogger) {
        toolLogger.info(localAssemblyContigs.count() + " chimeras indicating strand-switch-less breakpoints");

        final JavaPairRDD<ChimericAlignment, byte[]> chimeraAndSequence =
                localAssemblyContigs
                        .filter(tig ->
                                SimpleStrandSwitchVariantDetector.splitPairStrongEnoughEvidenceForCA(tig.alignmentIntervals.get(0),
                                        tig.alignmentIntervals.get(1), SimpleStrandSwitchVariantDetector.MORE_RELAXED_ALIGNMENT_MIN_MQ,
                                        0))
                        .mapToPair(tig -> convertAlignmentIntervalsToChimericAlignment(tig, refSequenceDictionary.getValue())).cache();

        final JavaRDD<VariantContext> annotatedBNDs =
                chimeraAndSequence
                        .mapToPair(pair -> new Tuple2<>(new NovelAdjacencyReferenceLocations(pair._1, pair._2, refSequenceDictionary.getValue()), pair._1))
                        .groupByKey()
                        .mapToPair(noveltyAndEvidence -> inferBNDType(noveltyAndEvidence, broadcastReference.getValue(), refSequenceDictionary.getValue()))
                        .flatMap(noveltyTypeAndEvidence ->
                                AnnotatedVariantProducer
                                        .produceAnnotatedBNDmatesVcFromNovelAdjacency(noveltyTypeAndEvidence._1,
                                                noveltyTypeAndEvidence._2._1, noveltyTypeAndEvidence._2._2, broadcastReference, refSequenceDictionary, sampleId).iterator());

        SVVCFWriter.writeVCF(annotatedBNDs.collect(), vcfOutputFileName.replace(".vcf", "_transBND.vcf"),
                refSequenceDictionary.getValue(), toolLogger);
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
