package org.broadinstitute.hellbender.tools.spark.sv;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

/**
 * Represents a collection of relevant information on one locally assembled contig from short reads.
 * TODO: should we move AlignmentRegion class here or keep it as top level class?
 */
class LocalAssemblyContig implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * This contig was assembled based on FASTQ records collected by {@link FindBadGenomicKmersSpark} and numbered by it.
     */
    final long assemblyID;

    /**
     * This together with {@link #assemblyID} completely determines which contig we are talking about.
     */
    final String contigID;

    final String seq;

    /**
     * Semi-qual score of bases that's available from fermi-lite but not from SGA.
     * May be {@code null} if assembler doesn't provide this feature
     */
    final byte[] quals;

    /**
     * Read count for the contigs.
     * TODO: fermi-lite provides support for this, not SGA
     */
    final int readCountSupport;
    static final int READ_COUNT_SUPPORT_DEFAULT = -1;

    private List<AlignmentRegion> alignmentRecords;

    /**
     * Construct a contig with assembly id, contig/vertex id, bases, but without number of reads generating the contig, or qual.
     */
    LocalAssemblyContig(final long assemblyID, final String contigID, final String seq){
        this(assemblyID, contigID, seq, null, READ_COUNT_SUPPORT_DEFAULT, null);
    }

    /**
     * Construct a contig with assembly id, contig/vertex id, bases, number of reads generating the contig, and qual.
     * Currently fermi-lite supports the last two features, not SGA.
     */
    LocalAssemblyContig(final long assemblyID, final String contigID, final String seq,
                        final List<AlignmentRegion> alignmentRegions,
                        final int readCountSupport, final byte[] quals){
        this.assemblyID = assemblyID;
        this.contigID   = contigID;
        this.seq        = seq;
        this.alignmentRecords = alignmentRegions;
        this.quals = (quals == null) ? null : quals;
        this.readCountSupport = readCountSupport;
    }

    LocalAssemblyContig(final long assemblyID, final String contigID, final String seq,
                        final List<AlignmentRegion> alignmentRegions){
        this(assemblyID, contigID, seq, alignmentRegions, READ_COUNT_SUPPORT_DEFAULT, null);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        LocalAssemblyContig that = (LocalAssemblyContig) o;

        if (assemblyID != that.assemblyID) return false;
        if (!contigID.equals(that.contigID)) return false;
        if (seq.equals(that.seq)) return false;
        if (readCountSupport != that.readCountSupport) return false;

        return (quals!=null && that.quals!=null) ? Arrays.equals(quals, that.quals) : (quals==null && that.quals==null);
    }

    @Override
    public int hashCode() {
        int result = (int) (assemblyID ^ (assemblyID >>> 32));
        result = 31 * result + contigID.hashCode();
        result = 31 * result + seq.hashCode();
        result = 31 * result + readCountSupport;
        return result;
    }

    List<AlignmentRegion> getAlignmentRegions(){
        return alignmentRecords;
    }
}
