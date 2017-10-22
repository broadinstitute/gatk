package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.FilterLongReadAlignmentsSAMSpark;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Created by valentin on 10/11/17.
 */
public interface SVHaplotype  {

    String REF_HAPLOTYPE_NAME = "ref";
    String ALT_HAPLOTYPE_NAME = "alt";

    /**
     * Returns the cigar of this haplotype versus the reference
     * @return
     */
    List<AlignmentInterval> getReferenceAlignmentIntervals();

    String getName();

    int getLength();

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

    <T> List<List<AlignmentInterval>> align(final Iterable<T> input, Function<T, byte[]> basesOf);

    default <T> List<AlignedContig> alignContigs(final Iterable<AlignedContig> contigs) {
        final List<String> names = Utils.stream(contigs).map(c -> c.contigName).collect(Collectors.toList());
        final List<byte[]> bases = Utils.stream(contigs).map(c -> c.contigSequence).collect(Collectors.toList());
        final List<List<AlignmentInterval>> intervals = align(bases, Function.identity());
        final List<AlignedContig> result = new ArrayList<>(intervals.size());
        final Set<String> haplotypeName = Collections.singleton(getName());
        for (int i = 0; i < intervals.size(); i++) {
            final List<List<AlignmentInterval>> bestCombos = FilterLongReadAlignmentsSAMSpark.pickBestConfigurations(names.get(i), intervals.get(i), haplotypeName);
            final AlignedContig alignedContig = new AlignedContig(names.get(i), bases.get(i), bestCombos.get(0), bestCombos.size() > 1);
            result.add(alignedContig);
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
}
