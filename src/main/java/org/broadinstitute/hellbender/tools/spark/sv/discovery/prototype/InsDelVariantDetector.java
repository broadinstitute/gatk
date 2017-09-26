package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import scala.Tuple2;

import java.util.List;


final class InsDelVariantDetector implements VariantDetectorFromLocalAssemblyContigAlignments {

    @Override
    public void inferSvAndWriteVCF(final JavaRDD<AlignedContig> localAssemblyContigs, final String vcfOutputFileName,
                                   final Broadcast<ReferenceMultiSource> broadcastReference, final SAMSequenceDictionary refSequenceDictionary,
                                   final Logger toolLogger){

        final JavaPairRDD<byte[], List<ChimericAlignment>> chimericAlignments =
                localAssemblyContigs
                        .mapToPair(tig -> convertAlignmentIntervalToChimericAlignment(tig,
                                StructuralVariationDiscoveryArgumentCollection.
                                        DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH));

        // usual business as in DiscoverVariantsFromContigAlignmentsSAMSpark#discoverVariantsAndWriteVCF()
        final JavaRDD<VariantContext> annotatedVariants =
                chimericAlignments
                        .flatMapToPair(DiscoverVariantsFromContigAlignmentsSAMSpark::discoverNovelAdjacencyFromChimericAlignments)
                        .groupByKey()
                        .mapToPair(noveltyAndEvidence -> DiscoverVariantsFromContigAlignmentsSAMSpark.inferType(noveltyAndEvidence._1, noveltyAndEvidence._2))
                        .map(noveltyTypeAndEvidence -> DiscoverVariantsFromContigAlignmentsSAMSpark.annotateVariant(noveltyTypeAndEvidence._1,
                                noveltyTypeAndEvidence._2._1, null, noveltyTypeAndEvidence._2._2, broadcastReference));

        SVVCFWriter.writeVCF(annotatedVariants.collect(), vcfOutputFileName, refSequenceDictionary, toolLogger);
    }

    /**
     * Delegates to {@link ChimericAlignment#parseOneContig(AlignedContig, int, boolean, boolean)}.
     */
    private static Tuple2<byte[], List<ChimericAlignment>> convertAlignmentIntervalToChimericAlignment (final AlignedContig contig,
                                                                                                        final int minAlignmentBlockSize) {

        return new Tuple2<>(contig.contigSequence,
                            ChimericAlignment.parseOneContig(contig, minAlignmentBlockSize, false, false));
    }
}
