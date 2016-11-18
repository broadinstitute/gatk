package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.SequenceUtil;
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

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.Strandedness.THREE_TO_FIVE;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.Strandedness.FIVE_TO_THREE;

/**
 * Internal caller for calling structural variants.
 */
class SVVariantCallerInternal implements Serializable {
    private static final long serialVersionUID = 1L;
    private static final Logger log = LogManager.getLogger(SVVariantCallerInternal.class);

    /**
     * First step in calling variants: parse all alignment records for a single locally-assembled contig and generate
     * chimeric alignments if available.
     *
     * @param input     made of ((assemblyId, contigId), ({alignmentRegions}, sequence)) of a signalling locally-assembled contig
     *
     * @return          the chimeric alignments of this sequence (empty if the sequence does not have any alignments)
     */
    static Iterable<ChimericAlignment> findChimericAlignments(Tuple2<Iterable<AlignmentRegion>, byte[]> input) {
        final byte[] contigSequence = input._2();
        final Iterable<AlignmentRegion> alignmentRegionsIterable = input._1();
        final List<AlignmentRegion> alignmentRegions = StreamSupport.stream(alignmentRegionsIterable.spliterator(), false).sorted(Comparator.comparing(a -> a.startInAssembledContig)).collect(Collectors.toList());
        return getChimericAlignmentsFromAlignmentRegions(contigSequence, alignmentRegions, SVConstants.DEFAULT_MIN_ALIGNMENT_LENGTH);
    }

    // TODO: this is looking at one chimeric alignment, i.e. a pair of alignment records at a time, potentially making calling complex variant difficult
    /**
     * Second step in calling variants: produce a single {@link BreakpointAllele} from a single {@link ChimericAlignment} of a contig.
     *
     * @param chimericAlignment chimeric alignment of a locally-assembled contig
     */
    static Tuple2<BreakpointAllele, ChimericAlignment> extractBreakpointAllele(final ChimericAlignment chimericAlignment) {
        return new Tuple2<>(new BreakpointAllele(chimericAlignment), chimericAlignment);
    }

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

        vcBuilder = vcBuilder.alleles(produceAlleles(broadcastReference.getValue(), contig, start, end)).id(produceVariantId(breakpointAllele));

        vcBuilder = updateAttributes(vcBuilder, start, end, assembledBreakpointsPerAllele);

        return vcBuilder.make();
    }





    @VisibleForTesting
    static List<ChimericAlignment> getChimericAlignmentsFromAlignmentRegions(final byte[] sequence, final List<AlignmentRegion> alignmentRegionList, final Integer minAlignLength) {
        if (alignmentRegionList.isEmpty()) {
            return new ArrayList<>();
        }
        final List<ChimericAlignment> results = new ArrayList<>(alignmentRegionList.size() - 1);
        final Iterator<AlignmentRegion> iterator = alignmentRegionList.iterator();
        final List<String> insertionAlignmentRegions = new ArrayList<>();
        if ( iterator.hasNext() ) {
            AlignmentRegion current = iterator.next();
            while (alignmentRegionMapQTooLow(current) && iterator.hasNext()) {
                current = iterator.next();
            }
            while ( iterator.hasNext() ) {
                final AlignmentRegion next = iterator.next();
                if (alignmentRegionIsTooSmall(current, next, minAlignLength)) {
                    continue;
                }

                if (nextAlignmentRegionInpairLowMapQ(current, next, minAlignLength)) {
                    if (iterator.hasNext()) {
                        insertionAlignmentRegions.add(next.toPackedString());
                        continue;
                    } else {
                        break;
                    }
                }

                final AlignmentRegion previous = current;
                current = next;

                final byte[] sequenceCopy = Arrays.copyOf(sequence, sequence.length);

                String homology = getHomology(current, previous, sequenceCopy);
                String insertedSequence = getInsertedSequence(current, previous, sequenceCopy);

                final ChimericAlignment chimericAlignment = new ChimericAlignment(previous, current, insertedSequence, homology, insertionAlignmentRegions);

                results.add(chimericAlignment);
            }
        }
        return results;
    }

    @VisibleForTesting
    static boolean alignmentRegionIsTooSmall(final AlignmentRegion current, final AlignmentRegion next, final Integer minAlignLength) {
        return current.referenceInterval.size() - SVVariantCallerUtils.overlapOnContig(current, next) < minAlignLength;
    }

    // TODO: think about how to test
    static boolean alignmentRegionMapQTooLow(final AlignmentRegion next) {
        return next.mapQual < 60;
    }

    @VisibleForTesting
    static boolean nextAlignmentRegionInpairLowMapQ(final AlignmentRegion current, final AlignmentRegion next, final Integer minAlignLength) {
        return alignmentRegionMapQTooLow(next) ||
                alignmentRegionIsTooSmall(next, current, minAlignLength) ||
                current.referenceInterval.contains(next.referenceInterval) ||
                next.referenceInterval.contains(current.referenceInterval);
    }

    @VisibleForTesting
    static String getHomology(final AlignmentRegion current, final AlignmentRegion previous, final byte[] sequenceCopy) {
        String homology = "";
        if (previous.endInAssembledContig >= current.startInAssembledContig) {
            final byte[] homologyBytes = Arrays.copyOfRange(sequenceCopy, current.startInAssembledContig - 1, previous.endInAssembledContig);
            if (previous.referenceInterval.getStart() > current.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(homologyBytes, 0, homologyBytes.length);
            }
            homology = new String(homologyBytes);
        }
        return homology;
    }

    @VisibleForTesting
    static String getInsertedSequence(final AlignmentRegion current, final AlignmentRegion previous, final byte[] sequenceCopy) {
        String insertedSequence = "";
        if (previous.endInAssembledContig < current.startInAssembledContig - 1) {

            final int insertionStart;
            final int insertionEnd;

            insertionStart = previous.endInAssembledContig + 1;
            insertionEnd = current.startInAssembledContig - 1;

            final byte[] insertedSequenceBytes = Arrays.copyOfRange(sequenceCopy, insertionStart - 1, insertionEnd);
            if (previous.referenceInterval.getStart() > current.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(insertedSequenceBytes, 0, insertedSequenceBytes.length);
            }
            insertedSequence = new String(insertedSequenceBytes);
        }
        return insertedSequence;
    }





    @VisibleForTesting
    static List<Allele> produceAlleles(final ReferenceMultiSource reference, final String contig, final int start, final int end) throws IOException {
        return new ArrayList<>(Arrays.asList(Allele.create(new String(reference.getReferenceBases(null, new SimpleInterval(contig, start, start)).getBases()), true), Allele.create("<INV>")));
    }

    @VisibleForTesting
    static String produceVariantId(final BreakpointAllele breakpointAllele) {
        // todo : hack for now for inversion only
        return (breakpointAllele.determineStrandedness()==FIVE_TO_THREE ? GATKSVVCFHeaderLines.INV_5_TO_3 : GATKSVVCFHeaderLines.INV_3_TO_5) + SVConstants.VARIANT_ID_FIELD_SEPARATER +
                breakpointAllele.leftJustifiedLeftBreakpoint.getContig() + SVConstants.VARIANT_ID_FIELD_SEPARATER +
                breakpointAllele.leftJustifiedLeftBreakpoint.getStart() + SVConstants.VARIANT_ID_FIELD_SEPARATER +
                breakpointAllele.leftJustifiedRightBreakpoint.getStart();
    }

    /**
     * Utility structs for extraction information from the consensus BreakpointAllele out of multiple ChimericAlignments,
     * to be later added to annotations of the VariantContext extracted.
     */
    static final class BreakpointAlleleAnnotations implements Serializable {
        private static final long serialVersionUID = 1L;

        final Integer minMQ;
        final Integer minAL;
        final String asmID;
        final String contigID;
        final List<String> insSeqMappings;

        BreakpointAlleleAnnotations(final ChimericAlignment chimericAlignment){
            minMQ = Math.min(chimericAlignment.regionWithLowerCoordOnContig.mapQual, chimericAlignment.regionWithHigherCoordOnContig.mapQual);
            minAL = Math.min(chimericAlignment.regionWithLowerCoordOnContig.referenceInterval.size(), chimericAlignment.regionWithHigherCoordOnContig.referenceInterval.size()) - SVVariantCallerUtils.overlapOnContig(chimericAlignment.regionWithLowerCoordOnContig, chimericAlignment.regionWithHigherCoordOnContig);;
            asmID = chimericAlignment.assemblyId;
            contigID = chimericAlignment.contigId;
            insSeqMappings = chimericAlignment.insertionMappings;
        }
    }

    @VisibleForTesting
    static VariantContextBuilder updateAttributes(VariantContextBuilder vcBuilder, final int start, final int end,
                                                  Tuple2<BreakpointAllele, Iterable<ChimericAlignment>> consensusAndAlignments) {

        // alignments should be sorted in a deterministic order for reproducibility
        final List<BreakpointAlleleAnnotations> annotations = StreamSupport.stream(consensusAndAlignments._2().spliterator(), false)
                .sorted((final ChimericAlignment o1, final ChimericAlignment o2) -> { // sort by assembly id, then sort by contig id
                    if (o1.assemblyId.equals(o2.assemblyId)) return o1.contigId.compareTo(o2.contigId);
                    else return o1.assemblyId.compareTo(o2.assemblyId);
                })
                .map(BreakpointAlleleAnnotations::new).collect(Collectors.toList());

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

        // todo: inversion specific attributes below, to be updated later in INSDEL calling
        vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.INV.toString());

        if (breakpointAllele.determineStrandedness() == FIVE_TO_THREE) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INV_5_TO_3, "");
        }

        if (breakpointAllele.determineStrandedness() == THREE_TO_FIVE) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INV_3_TO_5, "");
        }

        return vcBuilder;
    }
}
