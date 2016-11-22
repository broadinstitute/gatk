package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.Strandedness.SAME_STRAND;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.Strandedness.FIVE_TO_THREE;

/**
 * Internal caller for calling structural variants.
 */
class SVVariantConsensusCall implements Serializable {
    private static final long serialVersionUID = 1L;
    private static final Logger log = LogManager.getLogger(SVVariantConsensusCall.class);

    /**
     * Third step in calling variants: produce a VC from a {@link BreakpointAllele} (consensus among different assemblies if they all point to the same breakpoint).
     *
     * @param assembledBreakpointsPerAllele     consensus among different assemblies if they all point to the same breakpoint
     * @param broadcastReference                broadcasted reference
     * @throws IOException                      due to read operations on the reference
     */
    static VariantContext callVariantsFromConsensus(final Tuple2<BreakpointAllele, Iterable<ChimericAlignment>> assembledBreakpointsPerAllele,
                                                    final Broadcast<ReferenceMultiSource> broadcastReference) throws IOException {

        final BreakpointAllele breakpointAllele = assembledBreakpointsPerAllele._1;
        final String contig = breakpointAllele.leftJustifiedLeftBreakpoint.getContig();
        final int start = breakpointAllele.leftJustifiedLeftBreakpoint.getStart();
        final int end = breakpointAllele.leftJustifiedRightBreakpoint.getStart();

        VariantContextBuilder vcBuilder = new VariantContextBuilder().chr(contig).start(start).stop(end);

        vcBuilder = vcBuilder.alleles(produceAlleles(breakpointAllele, broadcastReference.getValue(), contig, start, end))
                             .id(produceVariantId(breakpointAllele));

        return updateAttributes(vcBuilder, start, end, assembledBreakpointsPerAllele).make();
    }


    @VisibleForTesting
    static List<Allele> produceAlleles(BreakpointAllele breakpointAllele, final ReferenceMultiSource reference, final String contig, final int start, final int end) throws IOException {
        final Allele refAllele = Allele.create(new String(reference.getReferenceBases(null, new SimpleInterval(contig, start, start)).getBases()), true);
        final Allele altAllele;
        if (breakpointAllele.determineStrandedness()!=SAME_STRAND) {
            altAllele = Allele.create(SVConstants.VCF_ALT_ALLELE_STRING_INV);
        } else {
            final boolean isNovelDisjoint = (breakpointAllele.leftJustifiedRightBreakpoint.getStart()-breakpointAllele.leftJustifiedLeftBreakpoint.getEnd()!=1); // reference intervals right next to each other
            final boolean isNovelAdjacency = !breakpointAllele.insertedSequence.isEmpty();
            if (isNovelAdjacency && isNovelDisjoint)
                altAllele = Allele.create(SVConstants.VCF_ALT_ALLELE_STRING_INDEL);
            else
                altAllele = isNovelDisjoint ? Allele.create(SVConstants.VCF_ALT_ALLELE_STRING_INS): Allele.create(SVConstants.VCF_ALT_ALLELE_STRING_DEL);
        }
        return new ArrayList<>(Arrays.asList(refAllele, altAllele));
    }

    @VisibleForTesting
    static String produceVariantId(final BreakpointAllele breakpointAllele) {

        // todo : hack for now for inversion only
        final BreakpointAllele.Strandedness strand = breakpointAllele.determineStrandedness();
        final String startString;
        if (strand==SAME_STRAND) {
            final boolean isNovelDisjoint = (breakpointAllele.leftJustifiedRightBreakpoint.getStart()-breakpointAllele.leftJustifiedLeftBreakpoint.getEnd()!=1); // reference intervals right next to each other
            final boolean isNovelAdjacency = !breakpointAllele.insertedSequence.isEmpty();
            if (isNovelAdjacency && isNovelDisjoint)
                startString = GATKSVVCFHeaderLines.SVTYPES.INSDEL.name();
            else
                startString = isNovelDisjoint ? GATKSVVCFHeaderLines.SVTYPES.INS.name(): GATKSVVCFHeaderLines.SVTYPES.DEL.name();
        } else {
            startString = (strand==FIVE_TO_THREE ? GATKSVVCFHeaderLines.INV_5_TO_3 : GATKSVVCFHeaderLines.INV_3_TO_5);
        }

        return  startString + SVConstants.VARIANT_ID_FIELD_SEPARATER +
                breakpointAllele.leftJustifiedLeftBreakpoint.getContig() + SVConstants.VARIANT_ID_FIELD_SEPARATER +
                breakpointAllele.leftJustifiedLeftBreakpoint.getStart() + SVConstants.VARIANT_ID_FIELD_SEPARATER +
                breakpointAllele.leftJustifiedRightBreakpoint.getStart();
    }

    /**
     * Utility structs for extraction information from the consensus BreakpointAllele out of multiple ChimericAlignments,
     * to be later added to annotations of the VariantContext extracted.
     */
    static final class BreakpointAlleleEvidenceAnnotations implements Serializable {
        private static final long serialVersionUID = 1L;

        final Integer minMQ;
        final Integer minAL;
        final String asmID;
        final String contigID;
        final List<String> insSeqMappings;

        BreakpointAlleleEvidenceAnnotations(final ChimericAlignment chimericAlignment){
            minMQ = Math.min(chimericAlignment.regionWithLowerCoordOnContig.mapQual, chimericAlignment.regionWithHigherCoordOnContig.mapQual);
            minAL = Math.min(chimericAlignment.regionWithLowerCoordOnContig.referenceInterval.size(), chimericAlignment.regionWithHigherCoordOnContig.referenceInterval.size())
                    - SVVariantCallerUtils.overlapOnContig(chimericAlignment.regionWithLowerCoordOnContig, chimericAlignment.regionWithHigherCoordOnContig);
            asmID = chimericAlignment.assemblyId;
            contigID = chimericAlignment.contigId;
            insSeqMappings = chimericAlignment.insertionMappings;
        }
    }

    @VisibleForTesting
    static VariantContextBuilder updateAttributes(VariantContextBuilder vcBuilder, final int start, final int end,
                                                  final Tuple2<BreakpointAllele, Iterable<ChimericAlignment>> consensusAndAlignments) {

        // alignments should be sorted in a deterministic order for reproducibility
        final List<BreakpointAlleleEvidenceAnnotations> annotations = StreamSupport.stream(consensusAndAlignments._2().spliterator(), false)
                .sorted((final ChimericAlignment o1, final ChimericAlignment o2) -> { // sort by assembly id, then sort by contig id
                    if (o1.assemblyId.equals(o2.assemblyId)) return o1.contigId.compareTo(o2.contigId);
                    else return o1.assemblyId.compareTo(o2.assemblyId);
                })
                .map(BreakpointAlleleEvidenceAnnotations::new).collect(Collectors.toList());

        vcBuilder = vcBuilder.attribute(VCFConstants.END_KEY, end)
                .attribute(GATKSVVCFHeaderLines.SVLEN, end - start)
                .attribute(GATKSVVCFHeaderLines.TOTAL_MAPPINGS, annotations.size())
                .attribute(GATKSVVCFHeaderLines.HQ_MAPPINGS, annotations.stream().filter(annotation -> annotation.minMQ == SVConstants.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD).count()) // todo: should use == or >=?
                .attribute(GATKSVVCFHeaderLines.MAPPING_QUALITIES, annotations.stream().map(annotation -> String.valueOf(annotation.minMQ)).collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)))
                .attribute(GATKSVVCFHeaderLines.ALIGN_LENGTHS, annotations.stream().map(annotation -> String.valueOf(annotation.minAL)).collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)))
                .attribute(GATKSVVCFHeaderLines.MAX_ALIGN_LENGTH, annotations.stream().map(annotation -> annotation.minAL).max(Comparator.naturalOrder()).orElse(0))
                .attribute(GATKSVVCFHeaderLines.ASSEMBLY_IDS, annotations.stream().map(annotation -> annotation.asmID).collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)))
                .attribute(GATKSVVCFHeaderLines.CONTIG_IDS, annotations.stream().map(annotation -> annotation.contigID).collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)));

        final BreakpointAllele breakpointAllele = consensusAndAlignments._1();

        if (!breakpointAllele.insertedSequence.isEmpty()) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, breakpointAllele.insertedSequence);
        }

        final List<String> insertionMappings = annotations.stream().map(annotation -> annotation.insSeqMappings).flatMap(List::stream).sorted().collect(Collectors.toList());

        if (!insertionMappings.isEmpty()) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE_MAPPINGS, insertionMappings.stream().collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)));
        }

        if (!breakpointAllele.homology.isEmpty()) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.HOMOLOGY, breakpointAllele.homology);
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.HOMOLOGY_LENGTH, breakpointAllele.homology.length());
        }

        // todo: sv type specific attributes below, to be generalized and refactored later
        final BreakpointAllele.Strandedness strand = breakpointAllele.determineStrandedness();
        final GATKSVVCFHeaderLines.SVTYPES type;
        if (strand==SAME_STRAND) {
            final boolean isNovelDisjoint = (breakpointAllele.leftJustifiedRightBreakpoint.getStart()-breakpointAllele.leftJustifiedLeftBreakpoint.getEnd()==1); // reference intervals right next to each other
            final boolean isNovelAdjacency = !breakpointAllele.insertedSequence.isEmpty();
            if (isNovelAdjacency && isNovelDisjoint)
                type = GATKSVVCFHeaderLines.SVTYPES.INSDEL;
            else
                type = isNovelDisjoint ? GATKSVVCFHeaderLines.SVTYPES.INS: GATKSVVCFHeaderLines.SVTYPES.DEL;

            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.SVTYPE, type.toString());
        } else {
            type = GATKSVVCFHeaderLines.SVTYPES.INV;
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.SVTYPE, type.toString());
            vcBuilder = vcBuilder.attribute( (strand == FIVE_TO_THREE) ? GATKSVVCFHeaderLines.INV_5_TO_3 : GATKSVVCFHeaderLines.INV_3_TO_5, "");
        }

        return vcBuilder;
    }
}
