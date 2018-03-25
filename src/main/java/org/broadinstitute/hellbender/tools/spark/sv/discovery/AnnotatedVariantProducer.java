package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.Logger;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleNovelAdjacencyAndChimericAlignmentEvidence;
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
     * @param breakendMates                     BND variants of mates to each other, assumed to be of size 2
     * @param simpleNovelAdjacencyAndChimericAlignmentEvidence novel adjacency and contig alignment evidence that induced the novel adjacency
     * @throws IOException                      due to reference retrieval
     */
    public static List<VariantContext> produceAnnotatedBNDmatesVcFromNovelAdjacency(final Tuple2<BreakEndVariantType, BreakEndVariantType> breakendMates,
                                                                                    final SimpleNovelAdjacencyAndChimericAlignmentEvidence simpleNovelAdjacencyAndChimericAlignmentEvidence,
                                                                                    final Broadcast<ReferenceMultiSource> broadcastReference,
                                                                                    final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                                                                    final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls,
                                                                                    final String sampleId)
            throws IOException {

        final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype = simpleNovelAdjacencyAndChimericAlignmentEvidence.getNovelAdjacencyReferenceLocations();
        final List<SimpleNovelAdjacencyAndChimericAlignmentEvidence.SimpleChimeraAndNCAMstring> contigEvidence = simpleNovelAdjacencyAndChimericAlignmentEvidence.getAlignmentEvidence();

        final VariantContext firstMate = produceAnnotatedVcFromInferredTypeAndRefLocations(novelAdjacencyAndAltHaplotype, breakendMates._1, contigEvidence,
                broadcastReference, broadcastSequenceDictionary, broadcastCNVCalls, sampleId);
        final VariantContext secondMate = produceAnnotatedVcFromInferredTypeAndRefLocations(novelAdjacencyAndAltHaplotype, breakendMates._2, contigEvidence,
                broadcastReference, broadcastSequenceDictionary, broadcastCNVCalls, sampleId);

        final VariantContextBuilder builder0 = new VariantContextBuilder(firstMate);
        builder0.attribute(GATKSVVCFConstants.BND_MATEID_STR, secondMate.getID());

        final VariantContextBuilder builder1 = new VariantContextBuilder(secondMate);
        builder1.attribute(GATKSVVCFConstants.BND_MATEID_STR, firstMate.getID());

        return Arrays.asList(builder0.make(), builder1.make());
    }

    /**
     * Produces a VC from a {@link NovelAdjacencyAndAltHaplotype}
     * (consensus among different assemblies if they all point to the same breakpoint).
     * @param inferredType                      inferred type of variant
     * @param broadcastReference                broadcast reference
     * @param broadcastSequenceDictionary       broadcast reference sequence dictionary
     * @param broadcastCNVCalls                 broadcast of external CNV calls (can be null)
     * @param sampleId                          sample identifier of the current sample
     * @throws IOException                      due to read operations on the reference
     */
    public static VariantContext produceAnnotatedVcFromInferredTypeAndRefLocations(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype,
                                                                                   final SvType inferredType,
                                                                                   final Iterable<SimpleNovelAdjacencyAndChimericAlignmentEvidence.SimpleChimeraAndNCAMstring> contigEvidence,
                                                                                   final Broadcast<ReferenceMultiSource> broadcastReference,
                                                                                   final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                                                                   final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls,
                                                                                   final String sampleId)
            throws IOException {

        final boolean variantIsBND = inferredType instanceof BreakEndVariantType;
        final int applicableStop;
        final SimpleInterval refLoc;
        if ( variantIsBND ) {
            final BreakEndVariantType breakEnd = (BreakEndVariantType) inferredType;
            if (breakEnd.isTheUpstreamMate()) {
                refLoc = novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc();
            } else {
                refLoc = novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc();
            }
            applicableStop = refLoc.getEnd(); // BND formatted variant shouldn't have END, this is just for having a valid "stop" for VC
        } else {
            refLoc = novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc();
            applicableStop = novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getEnd();
        }

        // basic information and attributes
        final VariantContextBuilder vcBuilder = new VariantContextBuilder()
                .chr(refLoc.getContig()).start(refLoc.getStart()).stop(applicableStop)
                .alleles(produceAlleles(refLoc, broadcastReference.getValue(), inferredType))
                .id(inferredType.getInternalVariantId())
                .attribute(GATKSVVCFConstants.SVTYPE, inferredType.toString());

        if (inferredType.getSVLength() != SVContext.NO_LENGTH)
            vcBuilder.attribute(GATKSVVCFConstants.SVLEN, inferredType.getSVLength());

        // attributes from complications
        inferredType.getTypeSpecificAttributes().forEach(vcBuilder::attribute);
        novelAdjacencyAndAltHaplotype.getComplication().toVariantAttributes().forEach(vcBuilder::attribute);

        // evidence used for producing the novel adjacency
        getEvidenceRelatedAnnotations(contigEvidence).forEach(vcBuilder::attribute);

        if ( ! variantIsBND )
            vcBuilder.attribute(VCFConstants.END_KEY, applicableStop);

        final byte[] altHaplotypeSequence = novelAdjacencyAndAltHaplotype.getAltHaplotypeSequence();
        if (altHaplotypeSequence != null && altHaplotypeSequence.length != 0)
            vcBuilder.attribute(GATKSVVCFConstants.SEQ_ALT_HAPLOTYPE, new String(altHaplotypeSequence));

        if (broadcastCNVCalls != null && (! variantIsBND) ) {
            final String cnvCallAnnotation = getExternalCNVCallAnnotation(refLoc, applicableStop, broadcastSequenceDictionary, broadcastCNVCalls, sampleId);
            if (! "".equals(cnvCallAnnotation)) {
                vcBuilder.attribute(GATKSVVCFConstants.EXTERNAL_CNV_CALLS, cnvCallAnnotation);
            }
        }

        return vcBuilder.make();
    }

    private static String getExternalCNVCallAnnotation(final SimpleInterval refLoc,
                                                       final int end,
                                                       final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                                       final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls,
                                                       final String sampleId) {
        final SVInterval variantInterval = new SVInterval(broadcastSequenceDictionary.getValue().getSequenceIndex(refLoc.getContig()),refLoc.getStart(), end);
        final SVIntervalTree<VariantContext> cnvCallTree = broadcastCNVCalls.getValue();
        final String cnvCallAnnotation =
                Utils.stream(cnvCallTree.overlappers(variantInterval))
                        .map(overlapper -> formatExternalCNVCallAnnotation(overlapper.getValue(), sampleId))
                        .collect(Collectors.joining(","));
        return cnvCallAnnotation;
    }

    private static String formatExternalCNVCallAnnotation(final VariantContext externalCNVCall, String sampleId) {
        return externalCNVCall.getID() + ":"
                + externalCNVCall.getGenotype(sampleId).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT) + ":"
                + externalCNVCall.getGenotype(sampleId).getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_QUALITY_FORMAT);
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
    public static List<Allele> produceAlleles(final SimpleInterval refLoc, final ReferenceMultiSource reference, final SvType svType)
            throws IOException {

        final byte[] refBases = reference.getReferenceBases(refLoc).getBases();

        return produceAlleles(refBases, svType);
    }

    @VisibleForTesting
    public static List<Allele> produceAlleles(final byte[] refBases, final SvType svType) {
        return new ArrayList<>(Arrays.asList(Allele.create(new String(refBases), true), svType.getAltAllele()));
    }

    public static List<VariantContext> annotateBreakpointBasedCallsWithImpreciseEvidenceLinks(final List<VariantContext> assemblyDiscoveredVariants,
                                                                                              final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                                                                              final ReadMetadata metadata,
                                                                                              final ReferenceMultiSource reference,
                                                                                              final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection parameters,
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

        ChimericContigAlignmentEvidenceAnnotations(final ChimericAlignment chimericAlignment, final String goodNonCanonicalMappingSATag){
            minMQ = Math.min(chimericAlignment.regionWithLowerCoordOnContig.mapQual,
                             chimericAlignment.regionWithHigherCoordOnContig.mapQual);
            minAL = Math.min(chimericAlignment.regionWithLowerCoordOnContig.referenceSpan.size(),
                             chimericAlignment.regionWithHigherCoordOnContig.referenceSpan.size())
                    - AlignmentInterval.overlapOnContig(chimericAlignment.regionWithLowerCoordOnContig,
                                                        chimericAlignment.regionWithHigherCoordOnContig);
            sourceContigName = chimericAlignment.sourceContigName;
            insSeqMappings = chimericAlignment.insertionMappings;
            this.goodNonCanonicalMappingSATag = goodNonCanonicalMappingSATag;
        }
    }

    @VisibleForTesting
    static Map<String, Object> getEvidenceRelatedAnnotations(final Iterable<SimpleNovelAdjacencyAndChimericAlignmentEvidence.SimpleChimeraAndNCAMstring> splitAlignmentEvidence) {

        final List<ChimericContigAlignmentEvidenceAnnotations> annotations =
                Utils.stream(splitAlignmentEvidence)
                        .sorted(Comparator.comparing(evidence -> evidence.simpleChimera.sourceContigName))
                        .map(evidence -> new ChimericContigAlignmentEvidenceAnnotations(evidence.simpleChimera, evidence.goodNonCanonicalMappingSATag))
                        .collect(Collectors.toList());

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

        final List<String> goodCanonicalMapping = annotations.stream().map(annotation -> annotation.goodNonCanonicalMappingSATag)
                .filter(s -> ! s.equals(AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME)).collect(Collectors.toList());
        if (!goodCanonicalMapping.isEmpty()) {
            attributeMap.put(GATKSVVCFConstants.CTG_GOOD_NONCANONICAL_MAPPING, goodCanonicalMapping.stream().collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
        }

        return attributeMap;
    }
}
