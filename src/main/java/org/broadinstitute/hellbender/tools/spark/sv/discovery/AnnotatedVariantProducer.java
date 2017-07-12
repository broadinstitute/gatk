package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.StringUtils;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Given identified pair of breakpoints for a simple SV and its supportive evidence, i.e. chimeric alignments,
 * produce an annotated {@link VariantContext}.
 */
class AnnotatedVariantProducer implements Serializable {
    private static final long serialVersionUID = 1L;

    // TODO: 12/12/16 does not handle translocation yet
    /**
     * Produces a VC from a {@link NovelAdjacencyReferenceLocations}
     * (consensus among different assemblies if they all point to the same breakpoint).
     * @param novelAdjacencyReferenceLocations  consensus among different assemblies if they all point to the same breakpoint
     * @param inferredType                      inferred type of variant
     * @param contigAlignments                  chimeric alignments from contigs used for generating this novel adjacency
     * @param broadcastReference                broadcasted reference
     *
     * @throws IOException                      due to read operations on the reference
     */
    static VariantContext produceAnnotatedVcFromNovelAdjacency(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations,
                                                               final SvType inferredType,
                                                               final Iterable<ChimericAlignment> contigAlignments,
                                                               final Broadcast<ReferenceMultiSource> broadcastReference)
            throws IOException {

        final String contig = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig();
        final int    start  = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd();
        final int    end    = novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();

        // basic information and attributes
        final VariantContextBuilder vcBuilder = new VariantContextBuilder()
                .chr(contig).start(start).stop(end)
                .alleles(produceAlleles(novelAdjacencyReferenceLocations, broadcastReference.getValue(), inferredType))
                .id(inferredType.getInternalVariantId())
                .attribute(VCFConstants.END_KEY, end)
                .attribute(GATKSVVCFHeaderLines.SVTYPE, inferredType.toString())
                .attribute(GATKSVVCFHeaderLines.SVLEN, inferredType.getSVLength());

        // attributes from complications
        inferredType.getTypeSpecificAttributes().forEach(vcBuilder::attribute);
        parseComplicationsAndMakeThemAttributeMap(novelAdjacencyReferenceLocations.complication).forEach(vcBuilder::attribute);

        // evidence used for producing the novel adjacency
        getEvidenceRelatedAnnotations(contigAlignments).forEach(vcBuilder::attribute);

        return vcBuilder.make();
    }

    // TODO: 12/13/16 again ignoring translocation
    @VisibleForTesting
    static List<Allele> produceAlleles(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations,
                                       final ReferenceMultiSource reference, final SvType SvType)
            throws IOException {

        final String contig = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig();
        final int start = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getStart();
        final byte[] refBases = reference.getReferenceBases(null, new SimpleInterval(contig, start, start)).getBases();

        return new ArrayList<>(Arrays.asList(Allele.create(new String(refBases), true), SvType.getAltAllele()));
    }

    /**
     * Not testing this because the complications are already tested in the NovelAdjacencyReferenceLocations class' own test,
     * more testing here would be actually testing VCBuilder.
     * @param breakpointComplications
     */
    private static Map<String, Object> parseComplicationsAndMakeThemAttributeMap(final BreakpointComplications breakpointComplications) {

        final Map<String, Object> attributeMap = new HashMap<>();

        if (!breakpointComplications.getInsertedSequenceForwardStrandRep().isEmpty()) {
            attributeMap.put(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, breakpointComplications.getInsertedSequenceForwardStrandRep());
        }

        if (!breakpointComplications.getHomologyForwardStrandRep().isEmpty()) {
            attributeMap.put(GATKSVVCFHeaderLines.HOMOLOGY, breakpointComplications.getHomologyForwardStrandRep());
            attributeMap.put(GATKSVVCFHeaderLines.HOMOLOGY_LENGTH, breakpointComplications.getHomologyForwardStrandRep().length());
        }

        if (breakpointComplications.hasDuplicationAnnotation()) {
            attributeMap.put(GATKSVVCFHeaderLines.DUP_REPEAT_UNIT_REF_SPAN, breakpointComplications.getDupSeqRepeatUnitRefSpan().toString());
            if(!breakpointComplications.getCigarStringsForDupSeqOnCtg().isEmpty()) {
                attributeMap.put(GATKSVVCFHeaderLines.DUP_SEQ_CIGARS,
                        StringUtils.join(breakpointComplications.getCigarStringsForDupSeqOnCtg(), VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
            }
            attributeMap.put(GATKSVVCFHeaderLines.DUPLICATION_NUMBERS,
                    new int[]{breakpointComplications.getDupSeqRepeatNumOnRef(), breakpointComplications.getDupSeqRepeatNumOnCtg()});
            if(breakpointComplications.isDupAnnotIsFromOptimization()) {
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
            minMQ = Math.min(chimericAlignment.regionWithLowerCoordOnContig.mapQual,
                    chimericAlignment.regionWithHigherCoordOnContig.mapQual);
            minAL = Math.min(chimericAlignment.regionWithLowerCoordOnContig.referenceSpan.size(),
                    chimericAlignment.regionWithHigherCoordOnContig.referenceSpan.size())
                    - SVVariantDiscoveryUtils.overlapOnContig(chimericAlignment.regionWithLowerCoordOnContig,
                    chimericAlignment.regionWithHigherCoordOnContig);
            sourceContigName = chimericAlignment.sourceContigName;
            insSeqMappings = chimericAlignment.insertionMappings;
        }
    }

    @VisibleForTesting
    static Map<String, Object> getEvidenceRelatedAnnotations(final Iterable<ChimericAlignment> splitAlignmentEvidence) {

        final List<NovelAdjacencyEvidenceAnnotations> annotations = Utils.stream(splitAlignmentEvidence)
                .sorted(Comparator.comparing(ca -> ca.sourceContigName))
                .map(NovelAdjacencyEvidenceAnnotations::new).collect(Collectors.toList());

        final Map<String, Object> attributeMap = new HashMap<>();
        attributeMap.put(GATKSVVCFHeaderLines.TOTAL_MAPPINGS,    annotations.size());
        attributeMap.put(GATKSVVCFHeaderLines.HQ_MAPPINGS,       annotations.stream().filter(annotation -> annotation.minMQ == StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD).count());// todo: should use == or >=?
        attributeMap.put(GATKSVVCFHeaderLines.MAPPING_QUALITIES, annotations.stream().map(annotation -> String.valueOf(annotation.minMQ)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        attributeMap.put(GATKSVVCFHeaderLines.ALIGN_LENGTHS,     annotations.stream().map(annotation -> String.valueOf(annotation.minAL)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        attributeMap.put(GATKSVVCFHeaderLines.MAX_ALIGN_LENGTH,  annotations.stream().map(annotation -> annotation.minAL).max(Comparator.naturalOrder()).orElse(0));
        attributeMap.put(GATKSVVCFHeaderLines.CONTIG_NAMES,      annotations.stream().map(annotation -> annotation.sourceContigName).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));

        final List<String> insertionMappings = annotations.stream().map(annotation -> annotation.insSeqMappings).flatMap(List::stream).sorted().collect(Collectors.toList());

        if (!insertionMappings.isEmpty()) {
            attributeMap.put(GATKSVVCFHeaderLines.INSERTED_SEQUENCE_MAPPINGS, insertionMappings.stream().collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        }
        return attributeMap;
    }
}
