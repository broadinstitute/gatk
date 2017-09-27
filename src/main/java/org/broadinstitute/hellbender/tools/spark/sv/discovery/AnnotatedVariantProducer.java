package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.StringUtils;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

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
                        novelAdjacencyReferenceLocations.complication, inferredType.get(0), null, contigAlignments, broadcastReference);

        final VariantContext secondMate =
                produceAnnotatedVcFromInferredTypeAndRefLocations(novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc, -1,
                        novelAdjacencyReferenceLocations.complication, inferredType.get(1), null, contigAlignments, broadcastReference);

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
     * @param altHaplotypeSeq                   alt haplotype sequence (could be null)
     * @param contigAlignments                  chimeric alignments from contigs used for generating this novel adjacency
     * @param broadcastReference                broadcasted reference
     *
     * @throws IOException                      due to read operations on the reference
     */
    static VariantContext produceAnnotatedVcFromInferredTypeAndRefLocations(final SimpleInterval refLoc, final int end,
                                                                            final BreakpointComplications breakpointComplications,
                                                                            final SvType inferredType,
                                                                            final byte[] altHaplotypeSeq,
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

        if (altHaplotypeSeq!=null)
            vcBuilder.attribute(GATKSVVCFConstants.SEQ_ALT_HAPLOTYPE, new String(altHaplotypeSeq));

        return vcBuilder.make();
    }

    public static VariantContext produceAnnotatedVcFromEvidenceTargetLink(final EvidenceTargetLink e,
                                                                          final SvType svType,
                                                                          final ReadMetadata metadata,
                                                                          final ReferenceMultiSource reference) {
        final String sequenceName = metadata.getContigName(e.getPairedStrandedIntervals().getLeft().getInterval().getContig());
        final int start = e.getPairedStrandedIntervals().getLeft().getInterval().midpoint();
        final int end = e.getPairedStrandedIntervals().getRight().getInterval().midpoint();
        try {
            final VariantContextBuilder builder = new VariantContextBuilder()
                    .chr(sequenceName)
                    .start(start)
                    .stop(end)
                    .id(svType.variantId)
                    .alleles(produceAlleles(new SimpleInterval(sequenceName, start, start), reference, svType))
                    .attribute(VCFConstants.END_KEY, end)
                    .attribute(GATKSVVCFConstants.SVLEN, svType.getSVLength())
                    .attribute(GATKSVVCFConstants.SVTYPE, svType.toString())
                    .attribute(GATKSVVCFConstants.IMPRECISE, true)
                    .attribute(GATKSVVCFConstants.CIPOS, produceCIInterval(start, e.getPairedStrandedIntervals().getLeft().getInterval()))
                    .attribute(GATKSVVCFConstants.CIEND, produceCIInterval(end, e.getPairedStrandedIntervals().getRight().getInterval()))
                    .attribute(GATKSVVCFConstants.READ_PAIR_SUPPORT, e.getReadPairs())
                    .attribute(GATKSVVCFConstants.SPLIT_READ_SUPPORT, e.getSplitReads());
            return builder.make();
        } catch (IOException e1) {
            throw new GATKException("error reading reference base for variant context " + svType.variantId, e1);
        }
    }

    /**
     * Produces the string representation of a VCF 4.2-style SV CI interval centered around 'point'.
     */
    @VisibleForTesting
    static String produceCIInterval(final int point, final SVInterval ciInterval) {
        Utils.validate(ciInterval.getStart() <= point && ciInterval.getEnd() >= point, "Interval must contain point");
        return String.join(",",
                String.valueOf(ciInterval.getStart() - point),
                String.valueOf(ciInterval.getEnd() - point));
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
            if (!breakpointComplications.getCigarStringsForDupSeqOnCtg().isEmpty()) {
                attributeMap.put(GATKSVVCFConstants.DUP_SEQ_CIGARS,
                        StringUtils.join(breakpointComplications.getCigarStringsForDupSeqOnCtg(), VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
            }
            attributeMap.put(GATKSVVCFConstants.DUPLICATION_NUMBERS,
                    new int[]{breakpointComplications.getDupSeqRepeatNumOnRef(), breakpointComplications.getDupSeqRepeatNumOnCtg()});
            if (breakpointComplications.isDupAnnotIsFromOptimization()) {
                attributeMap.put(GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, "");
            }

            if (breakpointComplications.getDupSeqStrandOnCtg() != null) {
                attributeMap.put(GATKSVVCFConstants.DUP_ORIENTATIONS,
                        breakpointComplications.getDupSeqStrandOnCtg().stream().map(Strand::toString).collect(Collectors.joining()));
            }
        }
        return attributeMap;
    }

    static VariantContext annotateWithImpreciseEvidenceLinks(final VariantContext variant,
                                                             final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                                             final SAMSequenceDictionary referenceSequenceDictionary,
                                                             final ReadMetadata metadata,
                                                             final int defaultUncertainty) {
        if (variant.getStructuralVariantType() == StructuralVariantType.DEL) {
            SVContext svc = SVContext.of(variant);
            final int padding = (metadata == null) ? defaultUncertainty : (metadata.getMaxMedianFragmentSize() / 2);
            PairedStrandedIntervals svcIntervals = svc.getPairedStrandedIntervals(metadata, referenceSequenceDictionary, padding);

            final Iterator<Tuple2<PairedStrandedIntervals, EvidenceTargetLink>> overlappers = evidenceTargetLinks.overlappers(svcIntervals);
            int readPairs = 0;
            int splitReads = 0;
            while (overlappers.hasNext()) {
                final Tuple2<PairedStrandedIntervals, EvidenceTargetLink> next = overlappers.next();
                readPairs += next._2.getReadPairs();
                splitReads += next._2.getSplitReads();
                overlappers.remove();
            }
            final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(variant);
            if (readPairs > 0) {
                variantContextBuilder.attribute(GATKSVVCFConstants.READ_PAIR_SUPPORT, readPairs);
            }
            if (splitReads > 0) {
                variantContextBuilder.attribute(GATKSVVCFConstants.SPLIT_READ_SUPPORT, splitReads);
            }

            return variantContextBuilder.make();
        } else {
            return variant;
        }
    }

    /**
     * Utility structs for extraction information from the consensus NovelAdjacencyReferenceLocations out of multiple ChimericAlignments,
     * to be later added to annotations of the VariantContext extracted.
     */
    private static final class NovelAdjacencyEvidenceAnnotations implements Serializable {
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

        final List<NovelAdjacencyEvidenceAnnotations> annotations =
                Utils.stream(splitAlignmentEvidence)
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
