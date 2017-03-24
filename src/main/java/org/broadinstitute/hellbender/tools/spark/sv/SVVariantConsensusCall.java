package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.tools.spark.sv.NovelAdjacencyReferenceLocations.EndConnectionType.FIVE_TO_THREE;

/**
 * Given identified pair of breakpoints for a simple SV and its supportive evidence, i.e. chimeric alignments,
 * produce an VariantContext.
 */
class SVVariantConsensusCall implements Serializable {
    private static final long serialVersionUID = 1L;

    // TODO: 12/12/16 does not handle translocation yet
    /**
     * Third step in calling variants: produce a VC from a {@link NovelAdjacencyReferenceLocations} (consensus among different assemblies if they all point to the same breakpoint).
     *
     * @param breakpointPairAndItsEvidence      consensus among different assemblies if they all point to the same breakpoint
     * @param broadcastReference                broadcasted reference
     * @throws IOException                      due to read operations on the reference
     */
    static VariantContext callVariantsFromConsensus(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> breakpointPairAndItsEvidence,
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
        if (endConnectionType == FIVE_TO_THREE) { // no strand switch happening, so no inversion
            if (start==end) { // something is inserted
                final boolean hasNoDupSeq = novelAdjacencyReferenceLocations.complication.dupSeqForwardStrandRep.isEmpty();
                final boolean hasNoInsertedSeq = novelAdjacencyReferenceLocations.complication.insertedSequenceForwardStrandRep.isEmpty();
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
                final boolean hasNoDupSeq = novelAdjacencyReferenceLocations.complication.dupSeqForwardStrandRep.isEmpty();
                final boolean hasNoInsertedSeq = novelAdjacencyReferenceLocations.complication.insertedSequenceForwardStrandRep.isEmpty();
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
            final SvType.TYPES x = SvType.TYPES.valueOf(type.toString());
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

        if (!novelAdjacencyReferenceLocations.complication.insertedSequenceForwardStrandRep.isEmpty()) {
            attributeMap.put(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, novelAdjacencyReferenceLocations.complication.insertedSequenceForwardStrandRep);
        }

        if (!novelAdjacencyReferenceLocations.complication.homologyForwardStrandRep.isEmpty()) {
            attributeMap.put(GATKSVVCFHeaderLines.HOMOLOGY, novelAdjacencyReferenceLocations.complication.homologyForwardStrandRep);
            attributeMap.put(GATKSVVCFHeaderLines.HOMOLOGY_LENGTH, novelAdjacencyReferenceLocations.complication.homologyForwardStrandRep.length());
        }

        if (!novelAdjacencyReferenceLocations.complication.dupSeqForwardStrandRep.isEmpty()) {
            attributeMap.put(GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, novelAdjacencyReferenceLocations.complication.dupSeqForwardStrandRep);
            attributeMap.put(GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, new int[]{novelAdjacencyReferenceLocations.complication.dupSeqRepeatNumOnRef, novelAdjacencyReferenceLocations.complication.dupSeqRepeatNumOnCtg});
        }
        return attributeMap;
    }

    /**
     * Utility structs for extraction information from the consensus NovelAdjacencyReferenceLocations out of multiple ChimericAlignments,
     * to be later added to annotations of the VariantContext extracted.
     */
    static final class BreakpointEvidenceAnnotations implements Serializable {
        private static final long serialVersionUID = 1L;

        final Integer minMQ;
        final Integer minAL;
        final String asmID;
        final String contigID;
        final List<String> insSeqMappings;

        BreakpointEvidenceAnnotations(final ChimericAlignment chimericAlignment){
            minMQ = Math.min(chimericAlignment.regionWithLowerCoordOnContig.mapQual, chimericAlignment.regionWithHigherCoordOnContig.mapQual);
            minAL = Math.min(chimericAlignment.regionWithLowerCoordOnContig.referenceInterval.size(), chimericAlignment.regionWithHigherCoordOnContig.referenceInterval.size())
                    - SVVariantCallerUtils.overlapOnContig(chimericAlignment.regionWithLowerCoordOnContig, chimericAlignment.regionWithHigherCoordOnContig);
            asmID = chimericAlignment.regionWithLowerCoordOnContig.assemblyId;
            contigID = chimericAlignment.regionWithLowerCoordOnContig.contigId;
            insSeqMappings = chimericAlignment.insertionMappings;
        }
    }

    @VisibleForTesting
    static Map<String, Object> getEvidenceRelatedAnnotations(final Iterable<ChimericAlignment> splitAlignmentEvidence) {

        final List<BreakpointEvidenceAnnotations> annotations = StreamSupport.stream(splitAlignmentEvidence.spliterator(), false)
                .sorted((final ChimericAlignment o1, final ChimericAlignment o2) -> { // sort by assembly id, then sort by contig id
                    if (o1.regionWithLowerCoordOnContig.assemblyId.equals(o2.regionWithLowerCoordOnContig.assemblyId)) return o1.regionWithLowerCoordOnContig.contigId.compareTo(o2.regionWithLowerCoordOnContig.contigId);
                    else return o1.regionWithLowerCoordOnContig.assemblyId.compareTo(o2.regionWithLowerCoordOnContig.assemblyId);
                })
                .map(BreakpointEvidenceAnnotations::new).collect(Collectors.toList());

        final Map<String, Object> attributeMap = new HashMap<>();
        attributeMap.put(GATKSVVCFHeaderLines.TOTAL_MAPPINGS, annotations.size());
        attributeMap.put(GATKSVVCFHeaderLines.HQ_MAPPINGS, annotations.stream().filter(annotation -> annotation.minMQ == SVConstants.CallingStepConstants.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD).count());// todo: should use == or >=?
        attributeMap.put(GATKSVVCFHeaderLines.MAPPING_QUALITIES, annotations.stream().map(annotation -> String.valueOf(annotation.minMQ)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        attributeMap.put(GATKSVVCFHeaderLines.ALIGN_LENGTHS, annotations.stream().map(annotation -> String.valueOf(annotation.minAL)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        attributeMap.put(GATKSVVCFHeaderLines.MAX_ALIGN_LENGTH, annotations.stream().map(annotation -> annotation.minAL).max(Comparator.naturalOrder()).orElse(0));
        attributeMap.put(GATKSVVCFHeaderLines.ASSEMBLY_IDS, annotations.stream().map(annotation -> annotation.asmID).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        attributeMap.put(GATKSVVCFHeaderLines.CONTIG_IDS, annotations.stream().map(annotation -> annotation.contigID).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));

        final List<String> insertionMappings = annotations.stream().map(annotation -> annotation.insSeqMappings).flatMap(List::stream).sorted().collect(Collectors.toList());

        if (!insertionMappings.isEmpty()) {
            attributeMap.put(GATKSVVCFHeaderLines.INSERTED_SEQUENCE_MAPPINGS, insertionMappings.stream().collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        }
        return attributeMap;
    }
}
