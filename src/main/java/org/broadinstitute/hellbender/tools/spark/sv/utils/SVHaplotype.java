package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.collections4.IterableUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigAlignmentsConfigPicker;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Common class of all haplotypes.
 * <h3>Naming conventions</h3>
 *
 * <p>
 *     By convention reference haplotypes (those that coincide with the reference sequence) are named as
 *     <code>ref[:XXX]</code>, where <code>XXX</code> indicates what reference haplotype this is when there is more than one
 *     contigouous sequence. Notice that <code>XXX</code> could contain further sub-sections using the ":" separator.
 *     For example an alternative translocation from somewhere on <code>chr1</code> to somwhere else in <code>chr17</code> would
 *     result in two reference haplotypes that could be named as: <code>ref:chr1, ref:chr17</code>.
 *     In the case of a long deletion or inversion on the same chromosome we might just use the break-point positions or
 *     relative locations on the genome: <code>ref:chr1:121211, ref:chr1:150000112</code>, or <code>ref:chr1:upstream, ref:chr1:downstream</code>.
 * </p>
 * <p>
 *     In contrast, those haplotypes that represent an alternative haplotypes would be formatted as <code>alt[:IDX][:XXX]</code>
 *     where <code>IDX</code> makes reference to the haplotype's index in multi-allelic variants (0-based) or another way to discern them
 *     (e.g. <code>alt:-101, alt:-100</code> for two alternative deletions of 101 and 100bp respectively). As with
 *     reference haplotypes <code>XXX</code> would be use to designate separate sequences within the same alternative allele.
 * </p>
 * <p>
 *     Those haplotypes that are neither reference or alternative but rather some arbitrary homologous sequence can have any name.
 *     For example this would included assembled contigs.
 * </p>
 *
 * </p>
 *
 */
public interface SVHaplotype {

    /**
     * Standard name for the haplotype that represents the reference.
     */
    String REF_HAPLOTYPE_NAME = "ref";

    /**
     * Standard name for the haplotype that represents the alternative allele.
     */
    String ALT_HAPLOTYPE_NAME = "alt";

    /**
     * Checks whether this haplotype is the reference haplotype.
     * @return true iff so.
     */
    default boolean isReference() {
        return getName().startsWith(REF_HAPLOTYPE_NAME);
    }

    /**
     * Checks whether this haplotype is part of an alternative allele.
     * @return true iff so.
     */
    default boolean isAlternative() {
        return getName().startsWith(ALT_HAPLOTYPE_NAME);
    }

    /**
     * Check whether this is reference nor alternative.
     * @return never {@code null}.
     */
    default boolean isNeitherReferenceNorAlternative() {
        return !isReference() && !isAlternative();
    }

    /**
     * Returns list of alignment intervals of this haplotype vs the reference.
     */
    List<AlignmentInterval> getReferenceAlignmentIntervals();

    String getName();

    /**
     * Length of the haplotype.
     * @return
     */
    int getLength();

    default boolean isContig() {
        return !isReference() && !isAlternative();
    }

    SimpleInterval getVariantLocation();

    /**
     * Copies bases from the haplotype into an array.
     * @param offset
     * @param whereTo
     * @param destOffset
     * @param length
     * @throws IllegalArgumentException if {@code whereTo} is {@code null} or the indeces and length passed are not correct.
     */
    void copyBases(final int offset, final byte[] whereTo, final int destOffset, final int length);

    default byte[] getBases(final int from, final int length) {
        final byte[] result = new byte[length];
        copyBases(from, result, 0, length);
        return result;
    }

    default byte[] getBases() {
        return getBases(0, getLength());
    }

    default AlignedContig toAlignedContig() {
        return new AlignedContig(this, this.getName());
    }

    <T> List<List<AlignmentInterval>> align(final Iterable<T> input, Function<T, byte[]> basesOf);

    default <T> List<AlignedContig> alignContigs(final Iterable<AlignedContig> contigs) {
        final List<String> names = Utils.stream(contigs).map(c -> c.contigName).collect(Collectors.toList());
        final List<byte[]> bases = Utils.stream(contigs).map(c -> c.contigSequence).collect(Collectors.toList());
        final List<List<AlignmentInterval>> intervals = align(bases, Function.identity());
        final List<AlignedContig> realignedContigs = IntStream.range(0, names.size())
                .mapToObj(i -> new AlignedContig(names.get(i), bases.get(i), intervals.get(i), false))
                .collect(Collectors.toList());
        final List<AlignedContig> result = new ArrayList<>(intervals.size());
        final Set<String> haplotypeName = Collections.singleton(getName());
        for (int i = 0; i < realignedContigs.size(); i++) {
            final List<AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings> mappings =
                     AssemblyContigAlignmentsConfigPicker.pickBestConfigurations(realignedContigs.get(i), haplotypeName, 0.0);
            if (mappings.isEmpty()) {
                result.add(realignedContigs.get(i));
            } else {
                final AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings topMapping = mappings.get(0);
                if (topMapping.getGoodMappings().isEmpty()) {
                    result.add(realignedContigs.get(i));
                } else if (topMapping.getGoodMappings().size() == realignedContigs.get(i).alignmentIntervals.size()) {
                    result.add(realignedContigs.get(i));
                } else {
                    result.add(new AlignedContig(realignedContigs.get(i).contigName, realignedContigs.get(i).contigSequence, topMapping.getGoodMappings(),
                            mappings.size() > 1));
                }
            }
        }
        return result;
    }

    int insertSize(final int start, final int end);


    /**
     * Returns the mapping reference span for this haplotype.
     * <p>
     *     It might return null indicating that the haplotype does not map to a single contig in the reference.
     * </p>
     *
     * @return {@code null} if the haplotype does not map to a single contig.
     */
    default SimpleInterval getReferenceSpan() {
        final List<AlignmentInterval> intervals = getReferenceAlignmentIntervals();
        if (intervals.isEmpty()) {
            return null;
        } else if (intervals.size() == 1) {
            return intervals.get(0).referenceSpan;
        } else {
            final SimpleInterval first = intervals.get(0).referenceSpan;
            int minStart = first.getStart();
            int maxEnd = first.getEnd();
            final String contig = first.getContig();
            for (int i = 1; i < intervals.size(); i++) {
                final SimpleInterval next = intervals.get(i).referenceSpan;
                if (next.getContig().equals(contig)) {
                    if (minStart > next.getStart()) {
                        minStart = next.getStart();
                    }
                    if (maxEnd < next.getEnd()) {
                        maxEnd = next.getEnd();
                    }
                } else {
                    return null;
                }
            }
            return new SimpleInterval(contig, minStart, maxEnd);
        }
    }

    default Cigar getCigar() {
        final List<AlignmentInterval> intervals = getReferenceAlignmentIntervals();
        if (intervals.isEmpty()) {
            return new Cigar();
        } else if (intervals.size() == 1) {
            return intervals.get(0).cigarAlongReference();
        } else {
            final List<AlignmentInterval> sortedByReferenceSpan = intervals.stream()
                    .sorted(Comparator.comparingInt(ai -> ai.referenceSpan.getStart()))
                    .collect(Collectors.toList());
            final AlignmentInterval first = intervals.get(0);
            final String contig = first.referenceSpan.getContig();
            final boolean forward = first.forwardStrand;
            final ArrayDeque<CigarElement> cigarElements = new ArrayDeque<>(intervals.size() * 3);
            for (final CigarElement element : first.cigarAlongReference()) {
                if (!element.getOperator().isClipping()) {
                    cigarElements.add(element);
                }
            }
            AlignmentInterval previous = first;
            for (int i = 1; i < sortedByReferenceSpan.size(); i++) {
                final AlignmentInterval next = sortedByReferenceSpan.get(i);
                if (next.forwardStrand == forward) {
                    final SimpleInterval nextRefInterval = next.referenceSpan;
                    if (nextRefInterval.getContig().equals(contig)) {
                        if (nextRefInterval.getStart() <= previous.referenceSpan.getEnd()) {
                            return new Cigar(); //TODO is it possible that both alignments overlapping bases would be exactly the same (compatible) I don't think that is true in practice.
                        } else if (forward && previous.endInAssembledContig >= next.startInAssembledContig) {
                            return new Cigar(); //TODO this could be fixed by reassiging oerlaipping bases in the contig to one of the two intervals behore hand, but this is not the place to fix that.
                        } else if (!forward && next.endInAssembledContig >= previous.startInAssembledContig) {
                            return new Cigar(); //TODO same as above.
                        } else {
                            final int deletionLength = next.referenceSpan.getStart() - previous.referenceSpan.getEnd() - 1;
                            if (deletionLength > 0) {
                                cigarElements.add(new CigarElement(deletionLength, CigarOperator.D));
                            }
                            final int insertionLength = forward ? next.startInAssembledContig - previous.endInAssembledContig - 1
                                                             : previous.startInAssembledContig - next.endInAssembledContig - 1;
                            if (insertionLength > 0) {
                                cigarElements.add(new CigarElement(insertionLength, CigarOperator.I));
                            }
                            //TODO It could be the case that we add an I and D operations in sequence here.... as we don't
                            //TODO know how to align these bases, I think this should be fixed somewherelese .
                            for (final CigarElement element : next.cigarAlongReference()) {
                                if (!element.getOperator().isClipping()) {
                                    cigarElements.add(element);
                                }
                            }
                        }
                    } else {
                        return new Cigar();
                    }
                } else {
                    return new Cigar();
                }
            }
            return new Cigar(new ArrayList<>(cigarElements));
        }
    }

    String getVariantId();
}
