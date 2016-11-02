package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.collections4.IterableUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;
import scala.Tuple3;
import scala.Tuple4;
import scala.Tuple5;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion.InversionType.INV_3_TO_5;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion.InversionType.INV_5_TO_3;

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
     * @param contigSequence                the assembled sequence
     * @param alignmentRegionsIterable      alignment records given by aligner for the sequence
     * @return                              the chimeric alignments of this sequence (empty if the sequence does not have any alignments)
     */
    static Iterable<ChimericAlignment> assembleBreakpointsFromAlignmentRegions(final byte[] contigSequence, final Iterable<AlignmentRegion> alignmentRegionsIterable) {
        final List<AlignmentRegion> alignmentRegions = IterableUtils.toList(alignmentRegionsIterable);
        if (alignmentRegions.size() > 1) { // todo: should remove this "if checking" and always sort when the input is guaranteed to be of size > 1
            alignmentRegions.sort(Comparator.comparing(a -> a.startInAssembledContig));
        }
        return getChimericAlignmentsFromAlignmentRegions(contigSequence, alignmentRegions, SVConstants.DEFAULT_MIN_ALIGNMENT_LENGTH);
    }

    // TODO: this is looking at one chimeric alignment, i.e. a pair of alignment records at a time, potentially making calling complex variant difficult
    /**
     * Second step in calling variants: produce a single {@link BreakpointAllele} from a single {@link ChimericAlignment} of a contig.
     *
     * @param breakpointIdAndAssembledBreakpoint    the inner Tuple2 is a pair of (asmID, contigID) for the contig that generated this chimeric alignment
     */
    static Tuple4<BreakpointAllele, String, String, ChimericAlignment> keyByBreakpointAllele(final Tuple2<Tuple2<String, String>, ChimericAlignment> breakpointIdAndAssembledBreakpoint) {
        final ChimericAlignment chimericAlignment = breakpointIdAndAssembledBreakpoint._2();
        // TODO: 11/4/16 INSDEL calling should add another kind of BreakpointAllele here
        final BreakpointAllele breakpointAllele = chimericAlignment.involvesStrandSwitch() ? new BreakpointAlleleInversion(chimericAlignment) : new BreakpointAllele(chimericAlignment);
        return new Tuple4<>(breakpointAllele, breakpointIdAndAssembledBreakpoint._1._1, breakpointIdAndAssembledBreakpoint._1._2, chimericAlignment);
    }

    /**
     * Third step in calling variants: produce a VC from a {@link BreakpointAllele} (consensus among different assemblies if they all point to the same breakpoint).
     *
     * @param assembledBreakpointsPerAllele     consensus among different assemblies if they all point to the same breakpoint
     * @param broadcastReference                broadcasted reference
     * @throws IOException                      due to read operations on the reference
     */
    static VariantContext getVariantFromBreakpointAlleleAlignments(final Tuple2<BreakpointAllele, Iterable<Tuple3<String, String, ChimericAlignment>>> assembledBreakpointsPerAllele,
                                                                   final Broadcast<ReferenceMultiSource> broadcastReference) throws IOException {

        final BreakpointAllele breakpointAllele = assembledBreakpointsPerAllele._1;
        final String contig = breakpointAllele.leftAlignedLeftBreakpoint.getContig();
        final int start = breakpointAllele.leftAlignedLeftBreakpoint.getStart();
        final int end = breakpointAllele.leftAlignedRightBreakpoint.getStart();
        VariantContextBuilder vcBuilder = new VariantContextBuilder().chr(contig).start(start).stop(end);

        vcBuilder = vcBuilder.alleles(produceAlleles(broadcastReference.getValue(), contig, start, end))
                             .id(produceVariantId(breakpointAllele));
        vcBuilder = updateAttributes(vcBuilder, breakpointAllele, assembledBreakpointsPerAllele._2, start, end);

        return vcBuilder.make();
    }

    // -----------------------------------------------------------------------------------------------
    // Group 1: from contig alignments to breakpoints
    // -----------------------------------------------------------------------------------------------

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
            while (treatAlignmentRegionAsInsertion(current) && iterator.hasNext()) {
                current = iterator.next();
            }
            while ( iterator.hasNext() ) {
                final AlignmentRegion next = iterator.next();
                if (alignmentRegionIsTooSmall(current, next, minAlignLength)) {
                    continue;
                }

                if (treatNextAlignmentRegionInPairAsInsertion(current, next, minAlignLength)) {
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
    static boolean treatAlignmentRegionAsInsertion(final AlignmentRegion next) {
        return next.mapQual < 60;
    }

    @VisibleForTesting
    static boolean treatNextAlignmentRegionInPairAsInsertion(final AlignmentRegion current, final AlignmentRegion next, final Integer minAlignLength) {
        return treatAlignmentRegionAsInsertion(next) ||
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

    // -----------------------------------------------------------------------------------------------
    // Group 2: from consensus (novel adjacency/disjoint) alleles to standard VC
    // -----------------------------------------------------------------------------------------------

    @VisibleForTesting
    static List<Allele> produceAlleles(final ReferenceMultiSource reference, final String contig, final int start, final int end) throws IOException {
        return new ArrayList<>(Arrays.asList(Allele.create(new String(reference.getReferenceBases(null, new SimpleInterval(contig, start, start)).getBases()), true), Allele.create("<INV>")));
    }

    @VisibleForTesting
    static String produceVariantId(final BreakpointAllele breakpointAllele) {
        final BreakpointAlleleInversion inversionAllele = (BreakpointAllele.BreakpointAlleleInversion) breakpointAllele;
        return inversionAllele.getInversionType().name() + SVConstants.VARIANT_ID_FIELD_SEPARATER +
                inversionAllele.leftAlignedLeftBreakpoint.getContig() + SVConstants.VARIANT_ID_FIELD_SEPARATER +
                inversionAllele.leftAlignedLeftBreakpoint.getStart() + SVConstants.VARIANT_ID_FIELD_SEPARATER +
                inversionAllele.leftAlignedRightBreakpoint.getStart();
    }

    /**
     * Auxiliary struct for extracting information from {@link BreakpointAllele.BreakpointAlleleInversion} to make the call.
     * todo: add documentation, especially which concept this aux is helping with
     */
    static final class AuxStruct implements Serializable {
        private static final long serialVersionUID = 1L;

        // 5 entries: (min mapping qualities of the two alignment regions that formed the chimeric alignment)
        //            (min alignment length of the two alignment regions that formed the chimeric alignment)
        //            (local assembly id)
        //            (contig id)
        //            (mapping information (collected from each chimeric alignment) of the inserted sequence)
        // reason for going with this way is such that the order will not be affected by how Spark distribute assemblies over executors
        final List<Tuple5<Integer, Integer, String, String, List<String>>> breakpointAlignmentInformation;
//        // following lists have the same length as input list to {@link collectInfo()}
//        final List<Integer> mqs;                // mapping qualities
//        final List<Integer> minAlignLengths;    // alignment length of the flanking regions of each identified breakpoint
//        final List<String> asmIds;              // local assembly ids
//        final List<String> assembledContigIds;  // contig ids

//        final List<String> insertionMappings;   // mapping information (collected from each chimeric alignment) of the inserted sequence

        final int numAssembledBreakpoints;            // todo: number of identified breakpoints for this contig
        final int highMqMappings;                     // number of breakpoints with high quality mappings, i.e. those with flanking mapping regions of a mapQ higher than a threshold
        final int maxMinAlignLength;                  // maximum alignment length of the flanking regions out of all identified breakpoints, i.e. max of {@code minAlignLengths}

        AuxStruct(final List<Tuple3<String, String, ChimericAlignment>> assembledBreakpoints) {

            breakpointAlignmentInformation = assembledBreakpoints.stream().map(entry ->{
                final ChimericAlignment chimericAlignment = entry._3();
                final int assembledBreakpointAlignmentLength =
                        Math.min(chimericAlignment.region1.referenceInterval.size(), chimericAlignment.region2.referenceInterval.size()) - SVVariantCallerUtils.overlapOnContig(chimericAlignment.region1, chimericAlignment.region2);
                return new Tuple5<>(Math.min(chimericAlignment.region1.mapQual, chimericAlignment.region2.mapQual), // mqs
                                    assembledBreakpointAlignmentLength,     // alignLengths
                                    entry._1(),                             // assemblyID
                                    entry._2(),                             // contigID
                                    chimericAlignment.insertionMappings);   // each chimeric alignment produces (a potentially empty) list of mappings for the inserted sequence
            }).sorted((final Tuple5<Integer, Integer, String, String, List<String>> o1, final Tuple5<Integer, Integer, String, String, List<String>> o2) -> { // sort by assembly id, then sort by contig id
                if (o1._3().equals(o2._3())) return o1._4().compareTo(o2._4());
                else return o1._3().compareTo(o2._3());
            }).collect(Collectors.toList());

            numAssembledBreakpoints = assembledBreakpoints.size();
            highMqMappings = (int) breakpointAlignmentInformation.stream().map(Tuple5::_1).filter(mq -> mq == SVConstants.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD).count();
            maxMinAlignLength = breakpointAlignmentInformation.stream().map(Tuple5::_2).max(Comparator.naturalOrder()).orElse(0);

//            final int numBreakpoints = assembledBreakpoints.size();
//            mqs = new ArrayList<>(numBreakpoints);
//            minAlignLengths = new ArrayList<>(numBreakpoints);
//            asmIds = new ArrayList<>(numBreakpoints);
//            assembledContigIds = new ArrayList<>(numBreakpoints);
//            insertionMappings = new ArrayList<>(numBreakpoints);
//            collectInfo(assembledBreakpoints);
//            // step over identified breakpoints from pair chimeric alignments
//            assembledBreakpoints.forEach(outterpair -> {
//                final ChimericAlignment chimericAlignment = outterpair._2;
//                mqs.add(Math.min(chimericAlignment.region1.mapQual, chimericAlignment.region2.mapQual));
//                minAlignLengths.add(Math.min(chimericAlignment.region1.referenceInterval.size(), chimericAlignment.region2.referenceInterval.size()) - SVVariantCallerUtils.overlapOnContig(chimericAlignment.region1, chimericAlignment.region2));
//                asmIds.add(outterpair._1._1);
//                assembledContigIds.add(chimericAlignment.contigId);
//                insertionMappings.addAll(chimericAlignment.insertionMappings);
//            });
//            // original "for loop" version
//            for (final Tuple2<Tuple2<String,String>, ChimericAlignment> assembledBreakpointPair : assembledBreakpoints) {
//                final ChimericAlignment chimericAlignment = assembledBreakpointPair._2;
//                numAssembledBreakpoints++;
//                final int assembledBreakpointMapq = Math.min(chimericAlignment.region1.mapQual, chimericAlignment.region2.mapQual);
//                mqs.add(assembledBreakpointMapq);
//                if (assembledBreakpointMapq == SVConstants.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD) {
//                    highMqMappings++;
//                }
//                final int assembledBreakpointAlignmentLength =
//                        Math.min(chimericAlignment.region1.referenceInterval.size(), chimericAlignment.region2.referenceInterval.size())
//                                - SVVariantCallerUtils.overlapOnContig(chimericAlignment.region1, chimericAlignment.region2);
//                minAlignLengths.add(assembledBreakpointAlignmentLength);
//                maxMinAlignLength = Math.max(maxMinAlignLength, assembledBreakpointAlignmentLength);
//                asmIds.add(assembledBreakpointPair._1._1);
//                assembledContigIds.add(chimericAlignment.contigId);
//                insertionMappings.addAll(chimericAlignment.insertionMappings);
//            }
        }
    }

    @VisibleForTesting
    static VariantContextBuilder updateAttributes(VariantContextBuilder vcBuilder, final BreakpointAllele breakpointAllele,
                                                  final Iterable<Tuple3<String, String, ChimericAlignment>> alignments,
                                                  final int start, final int end) {

        final BreakpointAlleleInversion inversionAllele = (BreakpointAlleleInversion) breakpointAllele;
        final AuxStruct auxStruct = new AuxStruct(IterableUtils.toList(alignments));

        vcBuilder = vcBuilder.attribute(VCFConstants.END_KEY, end)
                .attribute(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.INV.toString())
                .attribute(GATKSVVCFHeaderLines.SVLEN, end - start)
                .attribute(GATKSVVCFHeaderLines.TOTAL_MAPPINGS, auxStruct.numAssembledBreakpoints)
                .attribute(GATKSVVCFHeaderLines.HQ_MAPPINGS, auxStruct.highMqMappings)
                .attribute(GATKSVVCFHeaderLines.MAPPING_QUALITIES, auxStruct.breakpointAlignmentInformation.stream().map(Tuple5::_1).map(String::valueOf).collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)))
                .attribute(GATKSVVCFHeaderLines.ALIGN_LENGTHS, auxStruct.breakpointAlignmentInformation.stream().map(Tuple5::_2).map(String::valueOf).collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)))
                .attribute(GATKSVVCFHeaderLines.MAX_ALIGN_LENGTH, auxStruct.maxMinAlignLength)
                .attribute(GATKSVVCFHeaderLines.ASSEMBLY_IDS, auxStruct.breakpointAlignmentInformation.stream().map(Tuple5::_3).collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)))
                .attribute(GATKSVVCFHeaderLines.CONTIG_IDS, auxStruct.breakpointAlignmentInformation.stream().map(Tuple5::_4).map(s -> s.replace(" ", "_")).collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)));
//                .attribute(GATKSVVCFHeaderLines.MAPPING_QUALITIES, auxStruct.mqs.stream().map(String::valueOf).sorted().collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)))
//                .attribute(GATKSVVCFHeaderLines.ALIGN_LENGTHS, auxStruct.minAlignLengths.stream().map(String::valueOf).sorted().collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)))
//                .attribute(GATKSVVCFHeaderLines.ASSEMBLY_IDS, auxStruct.asmIds.stream().sorted().collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)))
//                .attribute(GATKSVVCFHeaderLines.CONTIG_IDS, auxStruct.assembledContigIds.stream().map(s -> s.replace(" ", "_")).sorted().collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)));


        if (inversionAllele.insertedSequence.length() > 0) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, inversionAllele.insertedSequence);
        }

        final List<String> insertionMappings = auxStruct.breakpointAlignmentInformation.stream().map(Tuple5::_5).flatMap(List::stream).collect(Collectors.toList());

        if (!insertionMappings.isEmpty()) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE_MAPPINGS, insertionMappings.stream().sorted().collect(Collectors.joining(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR)));
        }

        if (inversionAllele.homology.length() > 0) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.HOMOLOGY, inversionAllele.homology);
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.HOMOLOGY_LENGTH, inversionAllele.homology.length());
        }

        if (inversionAllele.getInversionType() == INV_5_TO_3) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INV_5_TO_3, "");
        }

        if (inversionAllele.getInversionType() == INV_3_TO_5) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INV_3_TO_5, "");
        }

        return vcBuilder;
    }
}
