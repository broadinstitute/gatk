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
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
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
public class AnnotatedVariantProducer implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Given novel adjacency and inferred BND variant types, produce annotated (and mate-connected) VCF BND records.
     * @param novelAdjacencyReferenceLocations  novel adjacency suggesting BND records
     * @param inferredType                      BND variants of mates to each other, assumed to be of size 2
     * @param contigAlignments                  chimeric alignments of supporting contig
     * @param broadcastReference                reference
     * @throws IOException
     */
    public static List<VariantContext> produceAnnotatedBNDmatesVcFromNovelAdjacency(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations,
                                                                                    final List<SvType> inferredType,
                                                                                    final Iterable<ChimericAlignment> contigAlignments,
                                                                                    final Broadcast<ReferenceMultiSource> broadcastReference)
            throws IOException {

        Utils.validateArg(inferredType.size() == 2,
                "Input novel adjacency doesn't seem to suggest mated BND records: \n" +
                        novelAdjacencyReferenceLocations.toString() + "\n" +
                        Utils.stream(contigAlignments).map(ChimericAlignment::onErrStringRep).collect(Collectors.toList()));

        final VariantContext firstMate =
                produceAnnotatedVcFromInferredTypeAndRefLocations(novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc, -1,
                        novelAdjacencyReferenceLocations.complication, inferredType.get(0), contigAlignments, broadcastReference);

        final VariantContext secondMate =
                produceAnnotatedVcFromInferredTypeAndRefLocations(novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc, -1,
                        novelAdjacencyReferenceLocations.complication, inferredType.get(1), contigAlignments, broadcastReference);

        final VariantContextBuilder builder0 = new VariantContextBuilder(firstMate);
        builder0.attribute(GATKSVVCFConstants.BND_MATEID_STR, secondMate.getID());

        final VariantContextBuilder builder1 = new VariantContextBuilder(secondMate);
        builder1.attribute(GATKSVVCFConstants.BND_MATEID_STR, firstMate.getID());

        return Arrays.asList(builder0.make(), builder1.make());
    }

    // TODO: 12/12/16 does not handle translocation yet
    /**
     * Produces a VC from a {@link NovelAdjacencyReferenceLocations}
     * (consensus among different assemblies if they all point to the same breakpoint).
     * @param refLoc                            corresponds to POS field of the returned VC, hence must be a point location.
     * @param end                               END of the VC, assumed to be < 0 if for BND formatted variant
     * @param breakpointComplications           complications associated with this breakpoint
     * @param inferredType                      inferred type of variant
     * @param contigAlignments                  chimeric alignments from contigs used for generating this novel adjacency
     * @param broadcastReference                broadcasted reference
     *
     * @throws IOException                      due to read operations on the reference
     */
    static VariantContext produceAnnotatedVcFromInferredTypeAndRefLocations(final SimpleInterval refLoc, final int end,
                                                                            final BreakpointComplications breakpointComplications,
                                                                            final SvType inferredType,
                                                                            final Iterable<ChimericAlignment> contigAlignments,
                                                                            final Broadcast<ReferenceMultiSource> broadcastReference)
            throws IOException {

        final int applicableEnd = end < 0 ? refLoc.getEnd() : end; // BND formatted variant shouldn't have END

        // basic information and attributes
        final VariantContextBuilder vcBuilder = new VariantContextBuilder()
                .chr(refLoc.getContig()).start(refLoc.getStart()).stop(applicableEnd)
                .alleles(produceAlleles(refLoc, broadcastReference.getValue(), inferredType))
                .id(inferredType.getInternalVariantId())
                .attribute(GATKSVVCFConstants.SVTYPE, inferredType.toString());

        if (inferredType.getSVLength() != SvType.INAPPLICABLE_LENGTH)
            vcBuilder.attribute(GATKSVVCFConstants.SVLEN, inferredType.getSVLength());

        // attributes from complications
        inferredType.getTypeSpecificAttributes().forEach(vcBuilder::attribute);
        parseComplicationsAndMakeThemAttributeMap(breakpointComplications).forEach(vcBuilder::attribute);

        // evidence used for producing the novel adjacency
        getEvidenceRelatedAnnotations(contigAlignments).forEach(vcBuilder::attribute);

        if (end > 0)
            vcBuilder.attribute(VCFConstants.END_KEY, applicableEnd);
        return vcBuilder.make();
    }

    // TODO: 12/13/16 again ignoring translocation
    @VisibleForTesting
    static List<Allele> produceAlleles(final SimpleInterval refLoc,
                                       final ReferenceMultiSource reference, final SvType SvType)
            throws IOException {

        final byte[] refBases = reference.getReferenceBases(null, refLoc).getBases();

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
            attributeMap.put(GATKSVVCFConstants.INSERTED_SEQUENCE, breakpointComplications.getInsertedSequenceForwardStrandRep());
        }

        if (!breakpointComplications.getHomologyForwardStrandRep().isEmpty()) {
            attributeMap.put(GATKSVVCFConstants.HOMOLOGY, breakpointComplications.getHomologyForwardStrandRep());
            attributeMap.put(GATKSVVCFConstants.HOMOLOGY_LENGTH, breakpointComplications.getHomologyForwardStrandRep().length());
        }

        if (breakpointComplications.hasDuplicationAnnotation()) {
            attributeMap.put(GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, breakpointComplications.getDupSeqRepeatUnitRefSpan().toString());
            if(!breakpointComplications.getCigarStringsForDupSeqOnCtg().isEmpty()) {
                attributeMap.put(GATKSVVCFConstants.DUP_SEQ_CIGARS,
                        StringUtils.join(breakpointComplications.getCigarStringsForDupSeqOnCtg(), VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
            }
            attributeMap.put(GATKSVVCFConstants.DUPLICATION_NUMBERS,
                    new int[]{breakpointComplications.getDupSeqRepeatNumOnRef(), breakpointComplications.getDupSeqRepeatNumOnCtg()});
            if(breakpointComplications.isDupAnnotIsFromOptimization()) {
                attributeMap.put(GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, "");
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
                    - AlignmentInterval.overlapOnContig(chimericAlignment.regionWithLowerCoordOnContig,
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
        attributeMap.put(GATKSVVCFConstants.TOTAL_MAPPINGS,    annotations.size());
        attributeMap.put(GATKSVVCFConstants.HQ_MAPPINGS,       annotations.stream().filter(annotation -> annotation.minMQ == StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD).count());// todo: should use == or >=?
        attributeMap.put(GATKSVVCFConstants.MAPPING_QUALITIES, annotations.stream().map(annotation -> String.valueOf(annotation.minMQ)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        attributeMap.put(GATKSVVCFConstants.ALIGN_LENGTHS,     annotations.stream().map(annotation -> String.valueOf(annotation.minAL)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        attributeMap.put(GATKSVVCFConstants.MAX_ALIGN_LENGTH,  annotations.stream().map(annotation -> annotation.minAL).max(Comparator.naturalOrder()).orElse(0));
        attributeMap.put(GATKSVVCFConstants.CONTIG_NAMES,      annotations.stream().map(annotation -> annotation.sourceContigName).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));

        final List<String> insertionMappings = annotations.stream().map(annotation -> annotation.insSeqMappings).flatMap(List::stream).sorted().collect(Collectors.toList());

        if (!insertionMappings.isEmpty()) {
            attributeMap.put(GATKSVVCFConstants.INSERTED_SEQUENCE_MAPPINGS, insertionMappings.stream().collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        }
        return attributeMap;
    }
}
