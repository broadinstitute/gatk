package org.broadinstitute.hellbender.tools.spark.sv;

import com.github.lindenb.jbwa.jni.AlnRgn;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.Objects;

class AlignmentRegion {

    final String contigId;
    final String assemblyId;
    final Cigar forwardStrandCigar;
    final boolean forwardStrand;
    final SimpleInterval referenceInterval;
    final int mapqual;
    final int startInAssembledContig;
    final int endInAssembledContig;
    final int assembledContigLength;
    final int mismatches;

    public AlignmentRegion(final String assemblyId, final String contigId, final AlnRgn alnRgn) {
        this.contigId = contigId;
        this.assemblyId = assemblyId;
        this.forwardStrand = alnRgn.getStrand() == '+';
        final Cigar alignmentCigar = TextCigarCodec.decode(alnRgn.getCigar());
        this.forwardStrandCigar = forwardStrand ? alignmentCigar : CigarUtils.invertCigar(alignmentCigar);
        this.referenceInterval = new SimpleInterval(alnRgn.getChrom(), (int) alnRgn.getPos() + 1, (int) (alnRgn.getPos() + forwardStrandCigar.getReferenceLength()));
        this.mapqual = alnRgn.getMQual();
        this.assembledContigLength = forwardStrandCigar.getReadLength() + getTotalHardClipping(forwardStrandCigar);
        this.startInAssembledContig = startOfAlignmentInContig(forwardStrandCigar);
        this.endInAssembledContig = endOfAlignmentInContig(assembledContigLength, forwardStrandCigar);
        this.mismatches = alnRgn.getNm();
    }

    public AlignmentRegion(final String assemblyId, final String contigId, final Cigar forwardStrandCigar, final boolean forwardStrand, final SimpleInterval referenceInterval, final int mapqual, final int startInAssembledContig, final int endInAssembledContig, final int mismatches) {
        this.contigId = contigId;
        this.assemblyId = assemblyId;
        this.forwardStrandCigar = forwardStrandCigar;
        this.forwardStrand = forwardStrand;
        this.referenceInterval = referenceInterval;
        this.mapqual = mapqual;
        this.startInAssembledContig = startInAssembledContig;
        this.endInAssembledContig = endInAssembledContig;
        this.assembledContigLength = forwardStrandCigar.getReadLength();
        this.mismatches = mismatches;
    }

    public AlignmentRegion(final GATKRead read) {
        this.assemblyId = null;
        this.contigId = read.getName();
        this.forwardStrand = ! read.isReverseStrand();
        this.forwardStrandCigar = forwardStrand ? read.getCigar() : CigarUtils.invertCigar(read.getCigar());
        this.referenceInterval = new SimpleInterval(read);
        this.assembledContigLength = forwardStrandCigar.getReadLength() + getTotalHardClipping(forwardStrandCigar);
        this.startInAssembledContig = startOfAlignmentInContig(forwardStrandCigar);
        this.endInAssembledContig = endOfAlignmentInContig(assembledContigLength, forwardStrandCigar);
        this.mapqual = read.getMappingQuality();
        if (read.hasAttribute("NM")) {
            this.mismatches = read.getAttributeAsInteger("NM");
        } else {
            this.mismatches = 0;
        }
    }

    public int overlapOnContig(final AlignmentRegion other) {
        return Math.max(0, Math.min(endInAssembledContig + 1, other.endInAssembledContig + 1) - Math.max(startInAssembledContig, other.startInAssembledContig));
    }

    private static int getTotalHardClipping(final Cigar cigar) {
        final List<CigarElement> cigarElements = cigar.getCigarElements();
        if (cigarElements.size() == 0) {
            return 0;
        }
        if (cigarElements.size() == 1) {
            return cigarElements.get(0).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(0).getLength() : 0;
        }
        return (cigarElements.get(0).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(0).getLength() : 0) +
                (cigarElements.get(cigarElements.size() - 1).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(cigarElements.size() - 1).getLength() : 0);
    }

    private static int startOfAlignmentInContig(final Cigar cigar) {
        return getNumClippedBases(true, cigar) + 1;
    }

    private static int endOfAlignmentInContig(final int assembledContigLength, final Cigar cigar) {
        return assembledContigLength - getNumClippedBases(false, cigar);
    }

    private static int getNumClippedBases(final boolean fromStart, final Cigar cigar) {
        int result = 0;
        int j = fromStart ? 0 : cigar.getCigarElements().size() - 1;
        final int offset = fromStart ? 1 : -1;
        CigarElement ce = cigar.getCigarElement(j);
        while (ce.getOperator().isClipping()) {
            result += ce.getLength();
            j += offset;
            if ( j < 0 || j >= cigar.getCigarElements().size() ) break;
            ce = cigar.getCigarElement(j);
        }
        return result;
    }

    @Override
    public String toString() {
        return assemblyId +
                "\t" +
                contigId +
                "\t" +
                referenceInterval.getContig() +
                "\t" +
                referenceInterval.getStart() +
                "\t" +
                referenceInterval.getEnd() +
                "\t" +
                (forwardStrand ? "+" : "-") +
                "\t" +
                forwardStrandCigar.toString() +
                "\t" +
                mapqual +
                "\t" +
                startInAssembledContig +
                "\t" +
                endInAssembledContig +
                "\t" +
                mismatches;
    }

    /**
     * Parses fields in the same order as they were output in toString:
     *
     * assemblyId
     * contigId
     * refContig
     * referenceInterval.getContig
     * referenceInterval.getStart
     * referenceInterval.getEnd
     * forwardStrand
     * forwardStrandCigar
     * startInAssembledContig
     * endInAssembledContig
     * mismatches
     *
     * @param fields
     * @return
     */
    public static AlignmentRegion fromString(final String[] fields) {
        final String breakpointId = fields[0];
        final String contigId = fields[1];
        final String refContig = fields[2];
        final Integer refStart = Integer.valueOf(fields[3]);
        final Integer refEnd = Integer.valueOf(fields[4]);
        final SimpleInterval refInterval = new SimpleInterval(refContig, refStart, refEnd);
        final boolean refStrand = ("+".equals(fields[5]));
        final Cigar cigar = TextCigarCodec.decode(fields[6]);
        final int mqual = Integer.valueOf(fields[7]);
        final int contigStart = Integer.valueOf(fields[8]);
        final int contigEnd = Integer.valueOf(fields[9]);
        final int mismatches = Integer.valueOf(fields[10]);
        return new AlignmentRegion(breakpointId, contigId, cigar, refStrand, refInterval, mqual, contigStart, contigEnd, mismatches);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        AlignmentRegion that = (AlignmentRegion) o;
        return forwardStrand == that.forwardStrand &&
                mapqual == that.mapqual &&
                startInAssembledContig == that.startInAssembledContig &&
                endInAssembledContig == that.endInAssembledContig &&
                assembledContigLength == that.assembledContigLength &&
                mismatches == that.mismatches &&
                Objects.equals(contigId, that.contigId) &&
                Objects.equals(assemblyId, that.assemblyId) &&
                Objects.equals(forwardStrandCigar, that.forwardStrandCigar) &&
                Objects.equals(referenceInterval, that.referenceInterval);
    }

    @Override
    public int hashCode() {
        return Objects.hash(contigId, assemblyId, forwardStrandCigar, forwardStrand, referenceInterval, mapqual, startInAssembledContig, endInAssembledContig, assembledContigLength, mismatches);
    }

    public String toPackedString() {
        return "" + startInAssembledContig + "-" + endInAssembledContig + ":" + referenceInterval.getContig() + ',' + referenceInterval.getStart() + ',' + (forwardStrand ? '+' : '-') + ',' + TextCigarCodec.encode(forwardStrandCigar) + ',' + mapqual + ',' + mismatches;
    }
}
