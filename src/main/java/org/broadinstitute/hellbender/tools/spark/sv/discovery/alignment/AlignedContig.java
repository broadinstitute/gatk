package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Locally assembled contig:
 *   its name
 *   its sequence as produced by the assembler (no reverse complement like in the SAM record if it maps to '-' strand), and
 *   its stripped-down alignment information.
 */
@DefaultSerializer(AlignedContig.Serializer.class)
public final class AlignedContig {

    private final String contigName;
    private final byte[] contigSequence;
    private final List<AlignmentInterval> alignmentIntervals;

    public AlignedContig(final GATKRead read) {
        Utils.nonNull(read);
        if (read.isUnmapped()) {
            throw new IllegalArgumentException("the input read cannot be unmapped");
        } else if (read.getCigar().isEmpty() || read.getCigar().containsOperator(CigarOperator.H)) {
            throw new IllegalArgumentException("the input read must have a cigar and cannot have hard-clips");
        } else {
            final byte[] bases = read.getBases();
            if (read.isReverseStrand()) {
                SequenceUtil.reverseComplement(bases);
            }
            final List<AlignmentInterval> intervals = new ArrayList<>();
            intervals.add(new AlignmentInterval(read));
            if (read.hasAttribute(SAMTag.SA.name())) {
                Arrays.stream(read.getAttributeAsString(SAMTag.SA.name()).split(";"))
                        .filter(s -> !s.isEmpty() && !s.equals("*"))
                        .map(AlignmentInterval::new)
                        .collect(Collectors.toList());
            }
            this.alignmentIntervals = Collections.unmodifiableList(intervals);
            this.contigSequence = bases;
            this.contigName = read.getName();
            this.hasEquallyGoodAlnConfigurations = false;
        }
    }

    /**
     * Composes a {@link AlignedContig} from all its constituted reads.
     * @param reads {@link Iterable} containing all relevant reads concerning the aligned-contig to compose.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code reads} is {@code null} or contains an inconsistent or incomplete
     *    set or read records.
     */
    public static AlignedContig of(final Iterable<GATKRead> reads) {
        Utils.nonNull(reads);
        return of(reads.iterator());
    }

    /**
     * Composes a {@link AlignedContig} from all its constituted reads.
     * @param reads {@link Stream} containing all the relevant read record for this contig.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code reads} is {@code null} or contains an inconsistent or incomplete
     *    set or read records.
     */
    public static AlignedContig of(final Stream<GATKRead> reads) {
        Utils.nonNull(reads);
        return of(reads.iterator());
    }

    /**
     * Composes a {@link AlignedContig} from all its constituted reads.
     * @param reads iterator that contains the reads as the remaining records to iterate upon.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code reads} is {@code null} or contains an inconsistent or incomplete
     *    set or read records.
     */
    public static AlignedContig of(final Iterator<GATKRead> reads) {
        Utils.nonNull(reads);
        if (!reads.hasNext()) {
            throw new IllegalArgumentException("there must be at least one read");
        } else {
            final GATKRead first = reads.next();
            if (!reads.hasNext()) {
                return new AlignedContig(first);
            } else {
                final String name = first.getName();
                final List<AlignmentInterval> alignmentIntervals = new ArrayList<>(5);
                final int length;
                final byte[] bases;
                if (first.isUnmapped()) {
                    length = first.getLength();
                    bases = Utils.nonNull(first.getBases(), "an unmapped record must have bases");
                    if (length != bases.length) {
                        throw new IllegalArgumentException("mismatch between bases array length and declared read length in unmapped record");
                    }
                } else {
                    final Cigar cigar = first.getCigar();
                    if (cigar == null || cigar.isEmpty()) {
                        throw new IllegalArgumentException("a mapped read must have a non-null nor empty cigar");
                    } else {
                        length = CigarUtils.countUnclippedReadBases(cigar);
                    }
                    bases = new byte[length];
                }
                GATKRead next = first;
                int basesDefined = 0;
                do {
                    if (next.isPaired()) {
                        throw new IllegalArgumentException("onlyt paired reads are currently supported");
                    } else if (!next.getName().equals(name)) {
                        throw new IllegalArgumentException("the input contains a mix of different reads");
                    } else {
                        basesDefined = mergeBases(basesDefined, bases, next);
                    }
                    if (!next.isUnmapped()) {
                        alignmentIntervals.add(new AlignmentInterval(next));
                    }
                    next = reads.hasNext() ? reads.next() : null;
                } while (next != null);
                if (basesDefined < length) {
                    throw new IllegalArgumentException("missing bases when looking across all the record provided");
                }
                return new AlignedContig(name, bases, alignmentIntervals, false);
            }
        }
    }

    private static int mergeBases(final int startBasesDefined, final byte[] bases, final GATKRead next) {
        if (startBasesDefined == bases.length) { // we have values for all the bases already we skip this step.
            return startBasesDefined;
        } else {
            final byte[] nextBases = Utils.nonNull(next.getBases(), "bases must be defined for every record");
            final int from, to;
            if (next.isUnmapped()) {
                if (nextBases.length != bases.length) {
                    throw new IllegalArgumentException("the input read bases for an unmapped read must be as long as the expected");
                }
                from = 0; to = bases.length;
            } else {
                final Cigar readCigar = next.getCigar();
                final Cigar cigar = next.isReverseStrand() ? CigarUtils.invertCigar(readCigar) : readCigar;
                if (cigar == null || cigar.isEmpty()) {
                    throw new IllegalArgumentException("mapped records must have a non-empty cigar");
                }
                from = CigarUtils.countLeftHardClippedBases(cigar);
                to = bases.length - CigarUtils.countRightHardClippedBases(cigar);
                if (next.isReverseStrand()) {
                    SequenceUtil.reverseComplement(nextBases);
                }
            }
            return startBasesDefined + mergeBases(nextBases, bases, from, to);
        }
    }

    private static int mergeBases(final byte[] src, final byte[] target, final int from, final int to) {
        int newlyDefined = 0;
        for (int i = from, j = 0; i < to; i++, j++) {
            final byte s = src[j], t = target[i];
            if (s == 0) {
                throw new IllegalArgumentException("a record contain null bases");
            } else if (t == 0) {
                target[i] = s;
                newlyDefined++;
            } else if (s != t) {
                throw new IllegalArgumentException("a record contains mismatching bases");
            }
        }
        return newlyDefined;
    }

    public AlignedContig(final Haplotype haplotype, final String name) {
        Utils.nonNull(haplotype);
        Utils.nonNull(name);
        final byte[] bases = haplotype.getBases();
        final Cigar cigar = haplotype.getCigar();
        Utils.validateArg(cigar != null && bases.length == cigar.getReadLength(), "invalid haplotype cigar");
        final Locatable loc = haplotype.getGenomeLocation();
        Utils.validateArg(loc != null && loc.getContig() != null && loc.getStart() <= loc.getEnd(), "invalid haplotype location");
        this.contigName = name;
        this.contigSequence = bases.clone();
        this.alignmentIntervals = Collections.singletonList(new AlignmentInterval(
                new SimpleInterval(loc.getContig(), loc.getStart(), loc.getStart() + cigar.getReferenceLength() - 1),
                        1, bases.length, cigar, false, SAMRecord.UNKNOWN_MAPPING_QUALITY,
                        AlignmentInterval.NO_NM, AlignmentInterval.NO_AS, false, false));
        this.hasEquallyGoodAlnConfigurations = false;
    }

    public AlignedContig(final String contigName, final byte[] contigSequence, final List<AlignmentInterval> alignmentIntervals,
                         final boolean hasEquallyGoodAlnConfigurations) {
        this.contigName = contigName;
        this.contigSequence = contigSequence;
        this.alignmentIntervals = alignmentIntervals.stream().sorted(getAlignmentIntervalComparator()).collect(Collectors.toList());
    }

    public AlignedContig(final Kryo kryo, final Input input) {

        contigName = input.readString();

        final int nBases = input.readInt();
        contigSequence = new byte[nBases];
        for (int b = 0; b < nBases; ++b) {
            contigSequence[b] = input.readByte();
        }

        final int nAlignments = input.readInt();
        alignmentIntervals = new ArrayList<>(nAlignments);
        for (int i = 0; i < nAlignments; ++i) {
            alignmentIntervals.add(new AlignmentInterval(kryo, input));
        }
    }

    static Comparator<AlignmentInterval> getAlignmentIntervalComparator() {
        Comparator<AlignmentInterval> comparePos = Comparator.comparingInt(aln -> aln.startInAssembledContig);
        Comparator<AlignmentInterval> compareRefTig = Comparator.comparing(aln -> aln.referenceSpan.getContig());
        Comparator<AlignmentInterval> compareRefSpanStart = Comparator.comparingInt(aln -> aln.referenceSpan.getStart());
        return comparePos.thenComparing(compareRefTig).thenComparing(compareRefSpanStart);
    }

    boolean hasOnly2Alignments() {
        return alignmentIntervals.size() == 2;
    }

    /**
     * @return first alignment of the contig
     */
    public AlignmentInterval getHeadAlignment() {
        return alignmentIntervals.get(0);
    }

    /**
     * @return last alignment of the contig
     */
    public AlignmentInterval getTailAlignment() {
        return alignmentIntervals.get(alignmentIntervals.size() - 1);
    }

    public String getContigName() {
        return contigName;
    }

    public byte[] getContigSequence() {
        return contigSequence;
    }

    public boolean isUnmapped() {
        return alignmentIntervals.isEmpty();
    }

    public List<AlignmentInterval> getAlignments() {
        return alignmentIntervals;
    }

    @Override
    public String toString() {
        return formatContigInfo(
                new Tuple2<>(contigName, alignmentIntervals.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList())));
    }

    /**
     * Format provided {@code tigNameAndMappings} for debugging.
     */
    public static String formatContigInfo(final Tuple2<String, List<String>> tigNameAndMappings) {
        return "(" + tigNameAndMappings._1 + ", " + tigNameAndMappings._2 + ")";
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final AlignedContig that = (AlignedContig) o;

        if (!contigName.equals(that.contigName)) return false;
        if (!Arrays.equals(contigSequence, that.contigSequence)) return false;
        return alignmentIntervals.equals(that.alignmentIntervals);
    }

    @Override
    public int hashCode() {
        int result = contigName.hashCode();
        result = 31 * result + Arrays.hashCode(contigSequence);
        result = 31 * result + alignmentIntervals.hashCode();
        return result;
    }

    void serialize(final Kryo kryo, final Output output) {

        output.writeString(contigName);

        output.writeInt(contigSequence.length);
        for (final byte base : contigSequence) {
            output.writeByte(base);
        }

        output.writeInt(alignmentIntervals.size());
        alignmentIntervals.forEach(it -> it.serialize(kryo, output));

    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AlignedContig> {
        @Override
        public void write(final Kryo kryo, final Output output, final AlignedContig alignedContig) {
            alignedContig.serialize(kryo, output);
        }

        @Override
        public AlignedContig read(final Kryo kryo, final Input input, final Class<AlignedContig> clazz) {
            return new AlignedContig(kryo, input);
        }
    }
}
