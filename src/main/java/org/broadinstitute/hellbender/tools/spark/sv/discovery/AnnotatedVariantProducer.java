package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.Logger;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleChimera;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleNovelAdjacencyAndChimericAlignmentEvidence;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection;

/**
 * Given identified pair of breakpoints for a simple SV and its supportive evidence, i.e. chimeric alignments,
 * produce an annotated {@link VariantContext}.
 */
public class AnnotatedVariantProducer implements Serializable {
    private static final long serialVersionUID = 1L;


    /**
     * Given novel adjacency and inferred variant types that should be linked together,
     * produce annotated, and linked VCF records.
     */
    public static List<VariantContext> produceLinkedAssemblyBasedVariants(final Tuple2<SvType, SvType> linkedVariants,
                                                                          final SimpleNovelAdjacencyAndChimericAlignmentEvidence simpleNovelAdjacencyAndChimericAlignmentEvidence,
                                                                          final Broadcast<ReferenceMultiSparkSource> broadcastReference,
                                                                          final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                                                          final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls,
                                                                          final String sampleId,
                                                                          final String linkKey) {

        final VariantContext firstVar = produceAnnotatedVcFromAssemblyEvidence(linkedVariants._1, simpleNovelAdjacencyAndChimericAlignmentEvidence,
                broadcastReference, broadcastSequenceDictionary, broadcastCNVCalls, sampleId).make();
        final VariantContext secondVar = produceAnnotatedVcFromAssemblyEvidence(linkedVariants._2, simpleNovelAdjacencyAndChimericAlignmentEvidence,
                broadcastReference, broadcastSequenceDictionary, broadcastCNVCalls, sampleId).make();

        final VariantContextBuilder builder1 = new VariantContextBuilder(firstVar);
        builder1.attribute(linkKey, secondVar.getID());

        final VariantContextBuilder builder2 = new VariantContextBuilder(secondVar);
        builder2.attribute(linkKey, firstVar.getID());

        // manually remove inserted sequence information from RPL event-produced DEL, when it can be linked with an INS
        if (linkedVariants._1 instanceof SimpleSVType.Deletion)
            return Arrays.asList(builder1.rmAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE)
                                         .rmAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE_LENGTH)
                                         .rmAttribute(GATKSVVCFConstants.SEQ_ALT_HAPLOTYPE)
                                         .make(),
                                 builder2.make());
        else if (linkedVariants._2 instanceof SimpleSVType.Deletion) {
            return Arrays.asList(builder1.make(),
                                 builder2.rmAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE)
                                         .rmAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE_LENGTH)
                                         .rmAttribute(GATKSVVCFConstants.SEQ_ALT_HAPLOTYPE)
                                         .make());
        } else
            return Arrays.asList(builder1.make(), builder2.make());
    }

    /**
     * Produces a VC from a {@link NovelAdjacencyAndAltHaplotype}
     * (consensus among different assemblies if they all point to the same breakpoint).
     * @param inferredType                      inferred type of variant
     * @param simpleNovelAdjacencyAndChimericAlignmentEvidence novel adjacency and evidence simple chimera contig(s) that support this {@code novelAdjacencyAndAltHaplotype}
     * @param broadcastReference                broadcast reference
     * @param broadcastSequenceDictionary       broadcast reference sequence dictionary
     * @param broadcastCNVCalls                 broadcast of external CNV calls, if available, will be used for annotating the assembly based calls
     * @param sampleId                          sample identifier of the current sample
     */
    public static VariantContextBuilder produceAnnotatedVcFromAssemblyEvidence(final SvType inferredType,
                                                                               final SimpleNovelAdjacencyAndChimericAlignmentEvidence simpleNovelAdjacencyAndChimericAlignmentEvidence,
                                                                               final Broadcast<ReferenceMultiSparkSource> broadcastReference,
                                                                               final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                                                               final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls,
                                                                               final String sampleId) {

        final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype = simpleNovelAdjacencyAndChimericAlignmentEvidence.getNovelAdjacencyReferenceLocations();
        final List<SimpleChimera> contigEvidence = simpleNovelAdjacencyAndChimericAlignmentEvidence.getAlignmentEvidence();

        // basic information and attributes
        final VariantContextBuilder vcBuilder = inferredType.getBasicInformation();

        // attributes from complications
        novelAdjacencyAndAltHaplotype.getComplication().toVariantAttributes().forEach(vcBuilder::attribute);

        // evidence used for producing the novel adjacency
        getAssemblyEvidenceRelatedAnnotations(contigEvidence).forEach(vcBuilder::attribute);

        // alt seq for non-BND variants, and if available or not empty
        final byte[] altHaplotypeSequence = novelAdjacencyAndAltHaplotype.getAltHaplotypeSequence();
        if (inferredType instanceof BreakEndVariantType) {
            return annotateWithExternalCNVCalls(inferredType.getVariantChromosome(), inferredType.getVariantStart(), inferredType.getVariantStop(),
                    vcBuilder, broadcastSequenceDictionary, broadcastCNVCalls, sampleId);
        } else if (altHaplotypeSequence != null && altHaplotypeSequence.length != 0)
            vcBuilder.attribute(GATKSVVCFConstants.SEQ_ALT_HAPLOTYPE, StringUtil.bytesToString(altHaplotypeSequence));

        return annotateWithExternalCNVCalls(inferredType.getVariantChromosome(), inferredType.getVariantStart(), inferredType.getVariantStop(),
                vcBuilder, broadcastSequenceDictionary, broadcastCNVCalls, sampleId);
    }

    public static VariantContext produceAnnotatedVcFromEvidenceTargetLink(final EvidenceTargetLink evidenceTargetLink,
                                                                          final SvType svType) {
        final PairedStrandedIntervals pairedStrandedIntervals = evidenceTargetLink.getPairedStrandedIntervals();
        final StrandedInterval strandedIntervalLeft = pairedStrandedIntervals.getLeft();
        final StrandedInterval strandedIntervalRight = pairedStrandedIntervals.getRight();
        final int start = strandedIntervalLeft.getInterval().midpoint();
        final int end = strandedIntervalRight.getInterval().midpoint();
        final VariantContextBuilder builder = svType
                .getBasicInformation()
                .attribute(GATKSVVCFConstants.CIPOS, produceCIInterval(start, strandedIntervalLeft.getInterval()))
                .attribute(GATKSVVCFConstants.CIEND, produceCIInterval(end, strandedIntervalRight.getInterval()))
                .attribute(GATKSVVCFConstants.READ_PAIR_SUPPORT, evidenceTargetLink.getReadPairs())
                .attribute(GATKSVVCFConstants.SPLIT_READ_SUPPORT, evidenceTargetLink.getSplitReads());
        return builder.make();
    }

    public static List<VariantContext> annotateBreakpointBasedCallsWithImpreciseEvidenceLinks(final List<VariantContext> assemblyDiscoveredVariants,
                                                                                              final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                                                                              final ReadMetadata metadata,
                                                                                              final ReferenceMultiSparkSource reference,
                                                                                              final DiscoverVariantsFromContigAlignmentsSparkArgumentCollection parameters,
                                                                                              final Logger localLogger) {

        final int originalEvidenceLinkSize = evidenceTargetLinks.size();
        final List<VariantContext> result = assemblyDiscoveredVariants
                .stream()
                .map(variant -> annotateWithImpreciseEvidenceLinks(
                        variant,
                        evidenceTargetLinks,
                        reference.getReferenceSequenceDictionary(null),
                        metadata, parameters.assemblyImpreciseEvidenceOverlapUncertainty))
                .collect(Collectors.toList());
        localLogger.info("Used " + (originalEvidenceLinkSize - evidenceTargetLinks.size()) + " evidence target links to annotate assembled breakpoints");
        return result;
    }

    //==================================================================================================================

    /**
     * Utility structs for extraction information from the consensus NovelAdjacencyAndAltHaplotype out of multiple ChimericAlignments,
     * to be later added to annotations of the VariantContext extracted.
     */
    private static final class ChimericContigAlignmentEvidenceAnnotations implements Serializable {
        private static final long serialVersionUID = 1L;

        final Integer minMQ;
        final Integer minAL;
        final String sourceContigName;
        final List<String> insSeqMappings;
        final String goodNonCanonicalMappingSATag;

        ChimericContigAlignmentEvidenceAnnotations(final SimpleChimera simpleChimera){
            minMQ = Math.min(simpleChimera.regionWithLowerCoordOnContig.mapQual,
                    simpleChimera.regionWithHigherCoordOnContig.mapQual);
            minAL = Math.min(simpleChimera.regionWithLowerCoordOnContig.referenceSpan.size(),
                    simpleChimera.regionWithHigherCoordOnContig.referenceSpan.size())
                    - AlignmentInterval.overlapOnContig(simpleChimera.regionWithLowerCoordOnContig,
                    simpleChimera.regionWithHigherCoordOnContig);
            sourceContigName = simpleChimera.sourceContigName;
            insSeqMappings = simpleChimera.insertionMappings;
            this.goodNonCanonicalMappingSATag = simpleChimera.goodNonCanonicalMappingSATag;
        }
    }

    @VisibleForTesting
    static Map<String, Object> getAssemblyEvidenceRelatedAnnotations(final Iterable<SimpleChimera> splitAlignmentEvidence) {

        final List<ChimericContigAlignmentEvidenceAnnotations> annotations =
                Utils.stream(splitAlignmentEvidence)
                        .sorted(Comparator.comparing(evidence -> evidence.sourceContigName))
                        .map(ChimericContigAlignmentEvidenceAnnotations::new)
                        .collect(Collectors.toList());

        final Map<String, Object> attributeMap = new HashMap<>();
        attributeMap.put(GATKSVVCFConstants.TOTAL_MAPPINGS,    annotations.size());
        attributeMap.put(GATKSVVCFConstants.HQ_MAPPINGS,       annotations.stream().filter(annotation -> annotation.minMQ == DiscoverVariantsFromContigAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD).count());
        attributeMap.put(GATKSVVCFConstants.MAPPING_QUALITIES, annotations.stream().map(annotation -> String.valueOf(annotation.minMQ)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        attributeMap.put(GATKSVVCFConstants.ALIGN_LENGTHS,     annotations.stream().map(annotation -> String.valueOf(annotation.minAL)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        attributeMap.put(GATKSVVCFConstants.MAX_ALIGN_LENGTH,  annotations.stream().map(annotation -> annotation.minAL).max(Comparator.naturalOrder()).orElse(0));
        attributeMap.put(GATKSVVCFConstants.CONTIG_NAMES,      annotations.stream().map(annotation -> annotation.sourceContigName).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));

        final List<String> insertionMappings = annotations.stream().map(annotation -> annotation.insSeqMappings).flatMap(List::stream).sorted().collect(Collectors.toList());

        if (!insertionMappings.isEmpty()) {
            attributeMap.put(GATKSVVCFConstants.INSERTED_SEQUENCE_MAPPINGS, insertionMappings.stream().collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        }

        final List<String> goodCanonicalMapping = annotations.stream().map(annotation -> annotation.goodNonCanonicalMappingSATag)
                .filter(s -> ! s.equals(AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME)).collect(Collectors.toList());
        if (!goodCanonicalMapping.isEmpty()) {
            attributeMap.put(GATKSVVCFConstants.CTG_GOOD_NONCANONICAL_MAPPING, goodCanonicalMapping.stream().collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        }

        return attributeMap;
    }

    @VisibleForTesting
    static VariantContextBuilder annotateWithExternalCNVCalls(final String recordContig, final int pos, final int end,
                                                              final VariantContextBuilder inputBuilder,
                                                              final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                                              final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls,
                                                              final String sampleId) {
        if (broadcastCNVCalls == null)
            return inputBuilder;
        final SVInterval variantInterval = new SVInterval(broadcastSequenceDictionary.getValue().getSequenceIndex(recordContig), pos, end);
        final SVIntervalTree<VariantContext> cnvCallTree = broadcastCNVCalls.getValue();
        final String cnvCallAnnotation =
                Utils.stream(cnvCallTree.overlappers(variantInterval))
                        .map(overlapper -> formatExternalCNVCallAnnotation(overlapper.getValue(), sampleId))
                        .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
        if (!cnvCallAnnotation.isEmpty()) {
            return inputBuilder.attribute(GATKSVVCFConstants.EXTERNAL_CNV_CALLS, cnvCallAnnotation);
        } else
            return inputBuilder;
    }

    private static String formatExternalCNVCallAnnotation(final VariantContext externalCNVCall, String sampleId) {
        return externalCNVCall.getID() + ":"
                + externalCNVCall.getGenotype(sampleId).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT) + ":"
                + externalCNVCall.getGenotype(sampleId).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_QUALITY_FORMAT);
    }

    //==================================================================================================================

    private static VariantContext annotateWithImpreciseEvidenceLinks(final VariantContext variant,
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
     * Produces the string representation of a VCF 4.2-style SV CI interval centered around 'point'.
     */
    @VisibleForTesting
    static String produceCIInterval(final int point, final SVInterval ciInterval) {
        Utils.validate(ciInterval.getStart() <= point && ciInterval.getEnd() >= point, "Interval must contain point");
        return String.join(",",
                String.valueOf(ciInterval.getStart() - point),
                String.valueOf(ciInterval.getEnd() - point));
    }
}
