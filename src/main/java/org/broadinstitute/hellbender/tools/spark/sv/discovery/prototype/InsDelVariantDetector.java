package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


final class InsDelVariantDetector implements VariantDetectorFromLocalAssemblyContigAlignments {

    @Override
    public void inferSvAndWriteVCF(final String vcfOutputFileName, final String sampleId, final JavaRDD<AlignedContig> contigs,
                                   final Broadcast<ReferenceMultiSource> broadcastReference,
                                   final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                   final Logger toolLogger){

        final JavaPairRDD<byte[], List<ChimericAlignment>> chimericAlignments =
                contigs
                        .mapToPair(tig -> convertAlignmentIntervalToChimericAlignment(tig,
                                StructuralVariationDiscoveryArgumentCollection.
                                        DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD,
                                StructuralVariationDiscoveryArgumentCollection.
                                        DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH,
                                broadcastSequenceDictionary.getValue()));

        // usual business as in DiscoverVariantsFromContigAlignmentsSAMSpark#discoverVariantsAndWriteVCF()
        final JavaRDD<VariantContext> annotatedVariants =
                chimericAlignments
                        .flatMapToPair(pair -> DiscoverVariantsFromContigAlignmentsSAMSpark.discoverNovelAdjacencyFromChimericAlignments(pair, broadcastSequenceDictionary.getValue()))
                        .groupByKey()
                        .mapToPair(noveltyAndEvidence -> DiscoverVariantsFromContigAlignmentsSAMSpark.inferType(noveltyAndEvidence._1, noveltyAndEvidence._2))
                        .map(noveltyTypeAndEvidence -> DiscoverVariantsFromContigAlignmentsSAMSpark.annotateVariant(noveltyTypeAndEvidence._1,
                                noveltyTypeAndEvidence._2._1, null, noveltyTypeAndEvidence._2._2, broadcastReference,
                                broadcastSequenceDictionary, null, sampleId));

        SVVCFWriter.writeVCF(annotatedVariants.collect(), vcfOutputFileName, broadcastSequenceDictionary.getValue(), toolLogger);
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
            if (!ca.isNotSimpleTranslocation())
                throw new GATKException.ShouldNeverReachHereException("Mapped assembled contigs are sent down the wrong path: " +
                        "contig suggesting \"translocation\" is sent down the insert/deletion path.\n" + contig.toString());
            results.add(ca);
        }
        return new Tuple2<>(contig.contigSequence,results);
    }
}
