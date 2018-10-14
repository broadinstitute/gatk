package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryInputMetaData;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;
import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH;

/**
 * This class scans the chimeric alignments of input {@link AlignedContig},
 * filters out the alignments that offers weak evidence for a breakpoint and,
 * makes interpretation based on the {@link SimpleChimera} extracted.
 */
public class ContigChimericAlignmentIterativeInterpreter {

    public static List<VariantContext> discoverVariantsFromChimeras(final SvDiscoveryInputMetaData svDiscoveryInputMetaData,
                                                                    final JavaRDD<AlignedContig> alignedContigs) {

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast =
                svDiscoveryInputMetaData.getReferenceData().getReferenceSequenceDictionaryBroadcast();

        // step 1: filter alignments and extract chimera pair
        final JavaPairRDD<byte[], List<SimpleChimera>> contigSeqAndChimeras =
                alignedContigs
                        .filter(alignedContig -> alignedContig.getAlignments().size() > 1)
                        .mapToPair(alignedContig -> {
                            final List<SimpleChimera> chimeras =
                                    parseOneContig(alignedContig, referenceSequenceDictionaryBroadcast.getValue(),
                                    true, DEFAULT_MIN_ALIGNMENT_LENGTH,
                                    CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, true);
                            return new Tuple2<>(alignedContig.getContigSequence(), chimeras);
                        });

        final Broadcast<ReferenceMultiSparkSource> referenceBroadcast = svDiscoveryInputMetaData.getReferenceData().getReferenceBroadcast();
        final List<SVInterval> assembledIntervals = svDiscoveryInputMetaData.getSampleSpecificData().getAssembledIntervals();
        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast = svDiscoveryInputMetaData.getSampleSpecificData().getCnvCallsBroadcast();
        final String sampleId = svDiscoveryInputMetaData.getSampleSpecificData().getSampleId();
        final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection discoverStageArgs = svDiscoveryInputMetaData.getDiscoverStageArgs();
        final Logger toolLogger = svDiscoveryInputMetaData.getToolLogger();

        // step 2: extract novel adjacency
        final JavaPairRDD<NovelAdjacencyAndAltHaplotype, Iterable<SimpleChimera>> narlsAndSources =
                contigSeqAndChimeras
                        .flatMapToPair(tigSeqAndChimeras -> {
                            final byte[] contigSeq = tigSeqAndChimeras._1;
                            final List<SimpleChimera> simpleChimeras = tigSeqAndChimeras._2;
                            final Stream<Tuple2<NovelAdjacencyAndAltHaplotype, SimpleChimera>> novelAdjacencyAndSourceChimera =
                                    simpleChimeras.stream()
                                            .map(ca -> new Tuple2<>(
                                                    new NovelAdjacencyAndAltHaplotype(ca, contigSeq,
                                                            referenceSequenceDictionaryBroadcast.getValue()), ca));
                            return novelAdjacencyAndSourceChimera.iterator();
                        })
                        .groupByKey()   // group the same novel adjacency produced by different contigs together
                        .cache();


        try {// step 3: evaluate performance turn into variant context
            SvDiscoveryUtils.evaluateIntervalsAndNarls(assembledIntervals, narlsAndSources.map(Tuple2::_1).collect(),
                    referenceSequenceDictionaryBroadcast.getValue(), discoverStageArgs, toolLogger);

            return narlsAndSources
                            .mapToPair(noveltyAndEvidence -> new Tuple2<>(inferSimpleTypeFromNovelAdjacency(noveltyAndEvidence._1, referenceBroadcast.getValue()),       // type inference based on novel adjacency and evidence alignments
                                    new SimpleNovelAdjacencyAndChimericAlignmentEvidence(noveltyAndEvidence._1, noveltyAndEvidence._2)))
                            .map(noveltyTypeAndEvidence ->
                                    AnnotatedVariantProducer
                                        .produceAnnotatedVcFromAssemblyEvidence(
                                                noveltyTypeAndEvidence._1, noveltyTypeAndEvidence._2,
                                                referenceBroadcast,
                                                referenceSequenceDictionaryBroadcast,
                                                cnvCallsBroadcast,
                                                sampleId).make()
                            )
                            .collect();
        } finally {
            narlsAndSources.unpersist();
        }
    }

    // =================================================================================================================

    /**
     * Parse all alignment records for a single locally-assembled contig and generate chimeric alignments if available.
     * Applies certain filters to skip the input alignment regions that are:
     *     1) if the alignment region's mapping quality is below a certain threshold, it is skipped
     *     2) if the alignment region is too small, it is skipped
     * If the alignment region passes the above two filters and the next alignment region could be treated as potential inserted sequence,
     * note down the mapping & alignment information of that region and skip it
     * @param alignedContig          made of (sorted {alignmentIntervals}, sequence) of a potentially-signalling locally-assembled contig
     * @param referenceDictionary    reference sequence dictionary
     * @param filterAlignmentByMqOrLength
     * @param uniqueRefSpanThreshold for an alignment interval to be used to construct a SimpleChimera,
     *                               how long a unique--i.e. the same ref span is not covered by other alignment intervals--alignment on the reference must it have
     * @param mapQualThresholdInclusive
     * @param filterWhollyContainedAlignments
     */
    @VisibleForTesting
    public static List<SimpleChimera> parseOneContig(final AlignedContig alignedContig,
                                                     final SAMSequenceDictionary referenceDictionary,
                                                     final boolean filterAlignmentByMqOrLength,
                                                     final int uniqueRefSpanThreshold,
                                                     final int mapQualThresholdInclusive,
                                                     final boolean filterWhollyContainedAlignments) {

        if (alignedContig.getAlignments().size() < 2) {
            return new ArrayList<>();
        }

        final Iterator<AlignmentInterval> iterator = alignedContig.getAlignments().iterator();

        // fast forward to the first alignment region with high MapQ
        AlignmentInterval current = iterator.next();
        if (filterAlignmentByMqOrLength) {
            while (mapQualTooLow(current, mapQualThresholdInclusive) && iterator.hasNext()) {
                current = iterator.next();
            }
        }

        final List<SimpleChimera> results = new ArrayList<>(alignedContig.getAlignments().size() - 1);
        final List<String> insertionMappings = new ArrayList<>();

        // then iterate over the AR's in pair to identify CA's.
        while ( iterator.hasNext() ) {
            final AlignmentInterval next = iterator.next();
            if (filterAlignmentByMqOrLength) {
                if (firstAlignmentIsTooShort(current, next, uniqueRefSpanThreshold)) {
                    continue;
                } else if (AssemblyContigAlignmentsConfigPicker.simpleChimeraWithStichableAlignments(current, next)) {
                    continue;
                } else if (nextAlignmentMayBeInsertion(current, next, mapQualThresholdInclusive, uniqueRefSpanThreshold, filterWhollyContainedAlignments)) {
                    if (iterator.hasNext()) {
                        insertionMappings.add(next.toPackedString());
                        continue;
                    } else {
                        break;
                    }
                }
            }

            // TODO: 10/18/17 this way of filtering CA based on not quality but alignment characteristics is temporary:
            //       this was initially developed for ins/del (and tested for that purpose), simple translocations travel through a different code path at the moment.
            // TODO: ultimately we need to merge these two code paths
            final SimpleChimera simpleChimera = new SimpleChimera(current, next, insertionMappings,
                    alignedContig.getContigName(), AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME,
                    referenceDictionary);
            // the following check/filter is due to the fact that simple translocations are to be handled in a different code path
            if (simpleChimera.isNeitherIncompleteNorSimpleTranslocation())
                results.add(simpleChimera);

            current = next;
        }

        return results;
    }

    // TODO: 11/22/16 it might also be suitable to consider the reference context this alignment region is mapped to
    //       and not simply apply a hard filter (need to think about how to test)
    private static boolean mapQualTooLow(final AlignmentInterval aln, final int mapQThresholdInclusive) {
        return aln.mapQual < mapQThresholdInclusive;
    }

    /**
     * @return if {@code first} is too short, when considering overlap with {@code second}
     */
    @VisibleForTesting
    public static boolean firstAlignmentIsTooShort(final AlignmentInterval first, final AlignmentInterval second,
                                                   final Integer minAlignLength) {
        return first.referenceSpan.size() - AlignmentInterval.overlapOnContig(first, second) < minAlignLength;
    }

    /**
     * To implement the idea that for two consecutive alignment regions of a contig, the one with higher reference coordinate might be a novel insertion.
     */
    @VisibleForTesting
    public static boolean nextAlignmentMayBeInsertion(final AlignmentInterval current, final AlignmentInterval next,
                                                      final Integer mapQThresholdInclusive, final Integer minAlignLength,
                                                      final boolean filterWhollyContained) {
        // not unique: inserted sequence may have low mapping quality (low reference uniqueness) or may be very small (low read uniqueness)
        final boolean isNotUnique = mapQualTooLow(next, mapQThresholdInclusive) || firstAlignmentIsTooShort(next, current, minAlignLength);
        return isNotUnique
                ||
                (filterWhollyContained && (current.referenceSpan.contains(next.referenceSpan) || next.referenceSpan.contains(current.referenceSpan)));
    }

    // =================================================================================================================

    @VisibleForTesting
    public static SimpleSVType inferSimpleTypeFromNovelAdjacency(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                                                                 final ReferenceMultiSparkSource reference) {

        final int start = novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getEnd();
        final int end = novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getStart();
        final StrandSwitch strandSwitch = novelAdjacencyAndAltHaplotype.getStrandSwitch();
        final boolean hasNoInsertedSeq = ! novelAdjacencyAndAltHaplotype.hasInsertedSequence();
        final boolean hasNoDupSeq = ! novelAdjacencyAndAltHaplotype.hasDuplicationAnnotation();

        final SimpleSVType type;
        if (strandSwitch == StrandSwitch.NO_SWITCH) { // no strand switch happening, so no inversion
            if (start==end) { // something is inserted
                if (hasNoDupSeq) {
                    if (hasNoInsertedSeq) {
                        throw new GATKException("Something went wrong in type inference, there's suspected insertion happening but no inserted sequence could be inferred "
                                + novelAdjacencyAndAltHaplotype.toString());
                    } else {
                        type = new SimpleSVType.Insertion(novelAdjacencyAndAltHaplotype, reference); // simple insertion (no duplication)
                    }
                } else {
                    type = new SimpleSVType.DuplicationTandem(novelAdjacencyAndAltHaplotype, reference);
                }
            } else {
                if (hasNoDupSeq) {
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.Deletion(novelAdjacencyAndAltHaplotype, reference); // clean deletion
                    } else {
                        type = new SimpleSVType.Deletion(novelAdjacencyAndAltHaplotype, reference); // scarred deletion
                    }
                } else {
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.Deletion(novelAdjacencyAndAltHaplotype, reference); // clean contraction of repeat 2 -> 1, or complex contraction
                    } else {
                        throw new GATKException("Something went wrong in novel adjacency interpretation: " +
                                " inferring simple SV type from a novel adjacency between two different reference locations, but annotated with both inserted sequence and duplication, which is NOT simple.\n"
                                + novelAdjacencyAndAltHaplotype.toString());
                    }
                }
            }
        } else {
            final int svLength = novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getStart() -
                    novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getEnd();
            type = new SimpleSVType.Inversion(novelAdjacencyAndAltHaplotype, svLength, reference);
        }

        return type;
    }
}
