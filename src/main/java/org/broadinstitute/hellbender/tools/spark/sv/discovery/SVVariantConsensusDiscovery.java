package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Given identified pair of breakpoints for a simple SV and its supportive evidence, i.e. chimeric alignments,
 * produce an VariantContext.
 */
public class SVVariantConsensusDiscovery implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Given contig alignments, scan for chimeric alignments which match a set of filtering criteria, and
     * emit novel adjacency not present on the reference used for aligning the contigs.
     */
    static JavaPairRDD<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> discoverNovelAdjacencyFromChimericAlignments(
            final JavaRDD<AlignedContig> alignedContigs,
            final Logger logger)
    {
        return alignedContigs.filter(alignedContig -> alignedContig.alignmentIntervals.size()>1) // filter out any contigs that has less than two alignment records
                .flatMapToPair(alignedContig -> NovelAdjacencyReferenceLocations.fromContigAlignments(alignedContig).iterator())
                .groupByKey();
    }

    // TODO: 12/12/16 does not handle translocation yet
    /**
     * Produces a VC from a {@link NovelAdjacencyReferenceLocations} (consensus among different assemblies if they all point to the same breakpoint).
     *
     * @param breakpointPairAndItsEvidence      consensus among different assemblies if they all point to the same breakpoint
     * @param broadcastReference                broadcasted reference
     * @throws IOException                      due to read operations on the reference
     */
    public static VariantContext discoverVariantsFromConsensus(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> breakpointPairAndItsEvidence,
                                                               final Broadcast<ReferenceMultiSource> broadcastReference)
            throws IOException {

        final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations = breakpointPairAndItsEvidence._1;
        final Iterable<ChimericAlignment> evidence = breakpointPairAndItsEvidence._2();
        final String contig = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig();
        final int start = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd();
        final int end = novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();

        Utils.validateArg(start<=end,
                "An identified breakpoint pair has left breakpoint positioned to the right of right breakpoint: " + novelAdjacencyReferenceLocations.toString());

        final SvType variant = getType(novelAdjacencyReferenceLocations);
        final VariantContextBuilder vcBuilder = new VariantContextBuilder()
                .chr(contig).start(start).stop(end)
                .alleles(produceAlleles(novelAdjacencyReferenceLocations, broadcastReference.getValue(), variant))
                .id(variant.getVariantId())
                .attribute(VCFConstants.END_KEY, end)
                .attribute(GATKSVVCFHeaderLines.SVTYPE, variant.toString())
                .attribute(GATKSVVCFHeaderLines.SVLEN, variant.getSVLength());

        variant.getTypeSpecificAttributes().forEach(vcBuilder::attribute);
        parseComplicationsAndMakeThemAttributeMap(novelAdjacencyReferenceLocations).forEach(vcBuilder::attribute);
        getEvidenceRelatedAnnotations(evidence).forEach(vcBuilder::attribute);

        return vcBuilder.make();
    }

    // TODO: 12/13/16 again ignoring translocation
    @VisibleForTesting
    static List<Allele> produceAlleles(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations, final ReferenceMultiSource reference, final SvType SvType)
            throws IOException {

        final String contig = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig();
        final int start = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getStart();

        final Allele refAllele = Allele.create(new String(reference.getReferenceBases(null, new SimpleInterval(contig, start, start)).getBases()), true);

        return new ArrayList<>(Arrays.asList(refAllele, SvType.getAltAllele()));
    }

    // TODO: 12/14/16 sv type specific attributes below, to be generalized and refactored later
    /**
     * Infer type of the variant, and produce an Id for the variant appropriately considering the complications.
     */
    @VisibleForTesting
    static SvType getType(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

        final int start = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd();
        final int end = novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
        final NovelAdjacencyReferenceLocations.EndConnectionType endConnectionType = novelAdjacencyReferenceLocations.endConnectionType;

        final SvType type;
        if (endConnectionType == NovelAdjacencyReferenceLocations.EndConnectionType.FIVE_TO_THREE) { // no strand switch happening, so no inversion
            if (start==end) { // something is inserted
                final boolean hasNoDupSeq = !novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation();
                final boolean hasNoInsertedSeq = novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep().isEmpty();
                if (hasNoDupSeq) {
                    if (hasNoInsertedSeq) {
                        throw new GATKException("Something went wrong in type inference, there's suspected insertion happening but no inserted sequence could be inferred " + novelAdjacencyReferenceLocations.toString());
                    } else {
                        type = new SvType.Insertion(novelAdjacencyReferenceLocations); // simple insertion (no duplication)
                    }
                } else {
                    if (hasNoInsertedSeq) {
                        type = new SvType.DuplicationTandem(novelAdjacencyReferenceLocations); // clean expansion of repeat 1 -> 2, or complex expansion
                    } else {
                        type = new SvType.DuplicationTandem(novelAdjacencyReferenceLocations); // expansion of 1 repeat on ref to 2 repeats on alt with inserted sequence in between the 2 repeats
                    }
                }
            } else {
                final boolean hasNoDupSeq = !novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation();
                final boolean hasNoInsertedSeq = novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep().isEmpty();
                if (hasNoDupSeq) {
                    if (hasNoInsertedSeq) {
                        type = new SvType.Deletion(novelAdjacencyReferenceLocations); // clean deletion
                    } else {
                        type = new SvType.Deletion(novelAdjacencyReferenceLocations); // scarred deletion
                    }
                } else {
                    if (hasNoInsertedSeq) {
                        type = new SvType.Deletion(novelAdjacencyReferenceLocations); // clean contraction of repeat 2 -> 1, or complex contraction
                    } else {
                        throw new GATKException("Something went wrong in type inference, there's suspected deletion happening but both inserted sequence and duplication exits (not supported yet): " + novelAdjacencyReferenceLocations.toString());
                    }
                }
            }
        } else {
            type = new SvType.Inversion(novelAdjacencyReferenceLocations);
        }

        // developer check to make sure new types are treated correctly
        try {
            SvType.TYPES.valueOf(type.toString());
        } catch (final IllegalArgumentException ex) {
            throw new GATKException.ShouldNeverReachHereException("Inferred type is not known yet: " + type.toString(), ex);
        }
        return type;
    }

    /**
     * Not testing this because the complications are already tested in the NovelAdjacencyReferenceLocations class' own test,
     * more testing here would be actually testing VCBuilder.
     */
    private static Map<String, Object> parseComplicationsAndMakeThemAttributeMap(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

        final Map<String, Object> attributeMap = new HashMap<>();

        if (!novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep().isEmpty()) {
            attributeMap.put(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep());
        }

        if (!novelAdjacencyReferenceLocations.complication.getHomologyForwardStrandRep().isEmpty()) {
            attributeMap.put(GATKSVVCFHeaderLines.HOMOLOGY, novelAdjacencyReferenceLocations.complication.getHomologyForwardStrandRep());
            attributeMap.put(GATKSVVCFHeaderLines.HOMOLOGY_LENGTH, novelAdjacencyReferenceLocations.complication.getHomologyForwardStrandRep().length());
        }

        if (novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation()) {
            attributeMap.put(GATKSVVCFHeaderLines.DUP_REPET_UNIT_REF_SPAN, novelAdjacencyReferenceLocations.complication.getDupSeqRepeatUnitRefSpan().toString());
            if(!novelAdjacencyReferenceLocations.complication.getCigarStringsForDupSeqOnCtg().isEmpty()) {
                attributeMap.put(GATKSVVCFHeaderLines.DUP_SEQ_CIGARS, StringUtils.join(novelAdjacencyReferenceLocations.complication.getCigarStringsForDupSeqOnCtg(), VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
            }
            attributeMap.put(GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, new int[]{novelAdjacencyReferenceLocations.complication.getDupSeqRepeatNumOnRef(), novelAdjacencyReferenceLocations.complication.getDupSeqRepeatNumOnCtg()});
            if(novelAdjacencyReferenceLocations.complication.isDupAnnotIsFromOptimization()) {
                attributeMap.put(GATKSVVCFHeaderLines.DUP_ANNOTATIONS_IMPRECISE, "");
            }
        }
        return attributeMap;
    }

    /**
     * Utility structs for extraction information from the consensus NovelAdjacencyReferenceLocations out of multiple ChimericAlignments,
     * to be later added to annotations of the VariantContext extracted.
     */
    static final class NovelAdjacencyEvidenceAnnotations implements Serializable {
        private static final long serialVersionUID = 1L;

        final Integer minMQ;
        final Integer minAL;
        final String sourceContigName;
        final List<String> insSeqMappings;

        NovelAdjacencyEvidenceAnnotations(final ChimericAlignment chimericAlignment){
            minMQ = Math.min(chimericAlignment.regionWithLowerCoordOnContig.mapQual, chimericAlignment.regionWithHigherCoordOnContig.mapQual);
            minAL = Math.min(chimericAlignment.regionWithLowerCoordOnContig.referenceInterval.size(), chimericAlignment.regionWithHigherCoordOnContig.referenceInterval.size())
                    - SVVariantDiscoveryUtils.overlapOnContig(chimericAlignment.regionWithLowerCoordOnContig, chimericAlignment.regionWithHigherCoordOnContig);
            sourceContigName = chimericAlignment.sourceContigName;
            insSeqMappings = chimericAlignment.insertionMappings;
        }
    }

    @VisibleForTesting
    static Map<String, Object> getEvidenceRelatedAnnotations(final Iterable<ChimericAlignment> splitAlignmentEvidence) {

        final List<NovelAdjacencyEvidenceAnnotations> annotations = StreamSupport.stream(splitAlignmentEvidence.spliterator(), false)
                .sorted(Comparator.comparing(ca -> ca.sourceContigName))
                .map(NovelAdjacencyEvidenceAnnotations::new).collect(Collectors.toList());

        final Map<String, Object> attributeMap = new HashMap<>();
        attributeMap.put(GATKSVVCFHeaderLines.TOTAL_MAPPINGS, annotations.size());
        attributeMap.put(GATKSVVCFHeaderLines.HQ_MAPPINGS, annotations.stream().filter(annotation -> annotation.minMQ == SVConstants.DiscoveryStepConstants.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD).count());// todo: should use == or >=?
        attributeMap.put(GATKSVVCFHeaderLines.MAPPING_QUALITIES, annotations.stream().map(annotation -> String.valueOf(annotation.minMQ)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        attributeMap.put(GATKSVVCFHeaderLines.ALIGN_LENGTHS, annotations.stream().map(annotation -> String.valueOf(annotation.minAL)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        attributeMap.put(GATKSVVCFHeaderLines.MAX_ALIGN_LENGTH, annotations.stream().map(annotation -> annotation.minAL).max(Comparator.naturalOrder()).orElse(0));
        attributeMap.put(GATKSVVCFHeaderLines.CONTIG_NAMES, annotations.stream().map(annotation -> annotation.sourceContigName).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));

        final List<String> insertionMappings = annotations.stream().map(annotation -> annotation.insSeqMappings).flatMap(List::stream).sorted().collect(Collectors.toList());

        if (!insertionMappings.isEmpty()) {
            attributeMap.put(GATKSVVCFHeaderLines.INSERTED_SEQUENCE_MAPPINGS, insertionMappings.stream().collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        }
        return attributeMap;
    }
}
