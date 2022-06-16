package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakEndVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SVInterval;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Simply a wrapper to link together
 * {@link NovelAdjacencyAndAltHaplotype} and
 * evidence {@link SimpleChimera}'s.
 */
@DefaultSerializer(SimpleNovelAdjacencyAndChimericAlignmentEvidence.Serializer.class)
public final class SimpleNovelAdjacencyAndChimericAlignmentEvidence {
    private static final NovelAdjacencyAndAltHaplotype.Serializer narlSerializer = new NovelAdjacencyAndAltHaplotype.Serializer();
    private static final SimpleChimera.Serializer alignmentEvidenceSerializer = new SimpleChimera.Serializer();

    private final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype;
    private final List<SimpleChimera> alignmentEvidence;

    public NovelAdjacencyAndAltHaplotype getNovelAdjacencyReferenceLocations() {
        return novelAdjacencyAndAltHaplotype;
    }
    public byte[] getAltHaplotypeSequence() {
        return novelAdjacencyAndAltHaplotype.getAltHaplotypeSequence();
    }
    public List<SimpleChimera> getAlignmentEvidence() {
        return alignmentEvidence;
    }

    public SimpleNovelAdjacencyAndChimericAlignmentEvidence(final NovelAdjacencyAndAltHaplotype novelAdjacencyReferenceLocations,
                                                            final Iterable<SimpleChimera> alignmentEvidence) {
        this.novelAdjacencyAndAltHaplotype = Utils.nonNull( novelAdjacencyReferenceLocations );
        this.alignmentEvidence = Lists.newArrayList( Utils.nonNull(alignmentEvidence) );
    }

    private SimpleNovelAdjacencyAndChimericAlignmentEvidence(final Kryo kryo, final Input input) {
        novelAdjacencyAndAltHaplotype = narlSerializer.read(kryo, input, NovelAdjacencyAndAltHaplotype.class);
        final int evidenceCount = input.readInt();
        alignmentEvidence = new ArrayList<>(evidenceCount);
        for (int i = 0; i < evidenceCount; ++i) {
            alignmentEvidence.add( alignmentEvidenceSerializer.read(kryo, input, SimpleChimera.class) );
        }
    }

    /**
     * This implementation is the 1st step going towards allowing re-interpretation,
     * below we simply take the inferred type and turn it to a VC,
     * future implementation may integrate other types of evidence and re-interpret if necessary
     */
    public List<VariantContext> toVariantContexts( final List<SvType> svTypes,
                                                   final String sampleId,
                                                   final SAMSequenceDictionary refDict,
                                                   final SVIntervalTree<VariantContext> cnvCalls) {
        if( svTypes.isEmpty() || svTypes.size() > 2 ) {
            throw new GATKException("Wrong number of variants sent for analysis: " + svTypes.toString() +
                    "\nWe currently only support 1 (symbolic simple or CPX) or 2 (BND mate pairs) variants for producing annotated variants.");
        }
        if (svTypes.size() == 2) {
            final SvType firstVar = svTypes.get(0);
            final SvType secondVar = svTypes.get(1);
            final String linkKey = firstVar instanceof BreakEndVariantType ? GATKSVVCFConstants.BND_MATEID_STR : GATKSVVCFConstants.LINK;
            return produceLinkedAssemblyBasedVariants(firstVar, secondVar, refDict, cnvCalls, sampleId, linkKey);
        } else {
            final VariantContext variantContext =
                    produceAnnotatedVcFromAssemblyEvidence(svTypes.get(0), refDict, cnvCalls, sampleId).make();
            return Collections.singletonList(variantContext);
        }
    }

    /**
     * Produces a VC from a {@link NovelAdjacencyAndAltHaplotype}
     * (consensus among different assemblies if they all point to the same breakpoint).
     * @param inferredType                      inferred type of variant
     * @param sequenceDictionary                reference sequence dictionary
     * @param cnvCalls                          external CNV calls, if available, will be used for annotating the assembly based calls
     * @param sampleId                          sample identifier of the current sample
     */
    @VisibleForTesting
    public VariantContextBuilder produceAnnotatedVcFromAssemblyEvidence( final SvType inferredType,
                                                                         final SAMSequenceDictionary sequenceDictionary,
                                                                         final SVIntervalTree<VariantContext> cnvCalls,
                                                                         final String sampleId) {

        // basic information and attributes
        final VariantContextBuilder vcBuilder = inferredType.getBasicInformation();

        // attributes from complications
        novelAdjacencyAndAltHaplotype.getComplication().toVariantAttributes().forEach(vcBuilder::attribute);

        // evidence used for producing the novel adjacency
        getAssemblyEvidenceRelatedAnnotations(alignmentEvidence).forEach(vcBuilder::attribute);

        // alt seq for non-BND variants, and if available or not empty
        final byte[] altHaplotypeSequence = novelAdjacencyAndAltHaplotype.getAltHaplotypeSequence();
        if (inferredType instanceof BreakEndVariantType) {
            return annotateWithExternalCNVCalls(inferredType.getVariantChromosome(), inferredType.getVariantStart(), inferredType.getVariantStop(),
                    vcBuilder, sequenceDictionary, cnvCalls, sampleId);
        }

        if (altHaplotypeSequence != null && altHaplotypeSequence.length != 0)
            vcBuilder.attribute(GATKSVVCFConstants.SEQ_ALT_HAPLOTYPE, StringUtil.bytesToString(altHaplotypeSequence));

        return annotateWithExternalCNVCalls(inferredType.getVariantChromosome(), inferredType.getVariantStart(), inferredType.getVariantStop(),
                vcBuilder, sequenceDictionary, cnvCalls, sampleId);
    }

    private static VariantContextBuilder annotateWithExternalCNVCalls(final String recordContig,
                                                                      final int pos,
                                                                      final int end,
                                                                      final VariantContextBuilder inputBuilder,
                                                                      final SAMSequenceDictionary sequenceDictionary,
                                                                      final SVIntervalTree<VariantContext> cnvCalls,
                                                                      final String sampleId) {
        if (cnvCalls == null)
            return inputBuilder;

        final SVInterval variantInterval = new SVInterval(sequenceDictionary.getSequenceIndex(recordContig), pos, end);
        final String cnvCallAnnotation =
                Utils.stream(cnvCalls.overlappers(variantInterval))
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

    /**
     * Given novel adjacency and inferred variant types that should be linked together,
     * produce annotated, and linked VCF records.
     */
    public List<VariantContext> produceLinkedAssemblyBasedVariants(final SvType svType1,
                                                                   final SvType svType2,
                                                                   final SAMSequenceDictionary sequenceDictionary,
                                                                   final SVIntervalTree<VariantContext> cnvCalls,
                                                                   final String sampleId,
                                                                   final String linkKey) {

        final VariantContext firstVar = produceAnnotatedVcFromAssemblyEvidence(
                    svType1, sequenceDictionary, cnvCalls, sampleId).make();
        final VariantContext secondVar = produceAnnotatedVcFromAssemblyEvidence(
                    svType2, sequenceDictionary, cnvCalls, sampleId).make();

        final VariantContextBuilder builder1 = new VariantContextBuilder(firstVar);
        builder1.attribute(linkKey, secondVar.getID());

        final VariantContextBuilder builder2 = new VariantContextBuilder(secondVar);
        builder2.attribute(linkKey, firstVar.getID());

        // manually remove inserted sequence information from RPL event-produced DEL, when it can be linked with an INS
        if (svType1 instanceof SimpleSVType.Deletion)
            return Arrays.asList(builder1.rmAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE)
                            .rmAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE_LENGTH)
                            .rmAttribute(GATKSVVCFConstants.SEQ_ALT_HAPLOTYPE)
                            .make(),
                    builder2.make());
        else if (svType2 instanceof SimpleSVType.Deletion) {
            return Arrays.asList(builder1.make(),
                    builder2.rmAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE)
                            .rmAttribute(GATKSVVCFConstants.INSERTED_SEQUENCE_LENGTH)
                            .rmAttribute(GATKSVVCFConstants.SEQ_ALT_HAPLOTYPE)
                            .make());
        } else
            return Arrays.asList(builder1.make(), builder2.make());
    }

    private static Map<String, Object> getAssemblyEvidenceRelatedAnnotations( final List<SimpleChimera> splitAlignmentEvidence ) {

        final List<ChimericContigAlignmentEvidenceAnnotations> annotations =
                splitAlignmentEvidence.stream()
                        .sorted(Comparator.comparing(evidence -> evidence.sourceContigName))
                        .map(ChimericContigAlignmentEvidenceAnnotations::new)
                        .collect(Collectors.toList());

        final Map<String, Object> attributeMap = new HashMap<>();
        attributeMap.put(GATKSVVCFConstants.TOTAL_MAPPINGS,    annotations.size());
        attributeMap.put(GATKSVVCFConstants.HQ_MAPPINGS,       annotations.stream().filter(annotation -> annotation.minMQ == StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD).count());
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

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SimpleNovelAdjacencyAndChimericAlignmentEvidence that = (SimpleNovelAdjacencyAndChimericAlignmentEvidence) o;

        if (!novelAdjacencyAndAltHaplotype.equals(that.novelAdjacencyAndAltHaplotype)) return false;
        return alignmentEvidence.equals(that.alignmentEvidence);
    }

    @Override
    public int hashCode() {
        int result = novelAdjacencyAndAltHaplotype.hashCode();
        result = 31 * result + alignmentEvidence.hashCode();
        return result;
    }

    private void serialize(final Kryo kryo, final Output output) {
        narlSerializer.write(kryo, output, novelAdjacencyAndAltHaplotype);
        output.writeInt(alignmentEvidence.size());
        alignmentEvidence.forEach(ev -> alignmentEvidenceSerializer.write(kryo, output, ev));
    }

    /**
     * Utility structs for extraction information from the consensus NovelAdjacencyAndAltHaplotype out of multiple ChimericAlignments,
     * to be later added to annotations of the VariantContext extracted.
     */
    public static final class ChimericContigAlignmentEvidenceAnnotations implements Serializable {
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

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SimpleNovelAdjacencyAndChimericAlignmentEvidence> {
        @Override
        public void write(final Kryo kryo, final Output output, final SimpleNovelAdjacencyAndChimericAlignmentEvidence novelAdjacencyReferenceLocations ) {
            novelAdjacencyReferenceLocations.serialize(kryo, output);
        }

        @Override
        public SimpleNovelAdjacencyAndChimericAlignmentEvidence read(final Kryo kryo, final Input input, final Class<SimpleNovelAdjacencyAndChimericAlignmentEvidence> klass ) {
            return new SimpleNovelAdjacencyAndChimericAlignmentEvidence(kryo, input);
        }
    }
}
