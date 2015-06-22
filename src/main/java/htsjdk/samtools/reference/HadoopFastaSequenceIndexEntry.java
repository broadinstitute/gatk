package htsjdk.samtools.reference;

/**
 * <b>Hellbender: Renamed copy of FastaSequenceIndexEntry. The only difference is that this class
 * is in its own source file to avoid Hellbender warnings. See https://github.com/broadinstitute/hellbender/issues/831.</b>
 * Hold an individual entry in a fasta sequence index file.
 */
class HadoopFastaSequenceIndexEntry {
    private String contig;
    private long location;
    private long size;
    private int basesPerLine;
    private int bytesPerLine;
    private final int sequenceIndex;

    /**
     * Create a new entry with the given parameters.
     * @param contig Contig this entry represents.
     * @param location Location (byte coordinate) in the fasta file.
     * @param size The number of bases in the contig.
     * @param basesPerLine How many bases are on each line.
     * @param bytesPerLine How many bytes are on each line (includes newline characters).
     */
    public HadoopFastaSequenceIndexEntry(String contig,
                                         long location,
                                         long size,
                                         int basesPerLine,
                                         int bytesPerLine,
                                         int sequenceIndex) {
        this.contig = contig;
        this.location = location;
        this.size = size;
        this.basesPerLine = basesPerLine;
        this.bytesPerLine = bytesPerLine;
        this.sequenceIndex = sequenceIndex;
    }

    /**
     * Gets the contig associated with this entry.
     * @return String representation of the contig.
     */
    public String getContig() {
        return contig;
    }

    /**
     * Sometimes contigs need to be adjusted on-the-fly to
     * match sequence dictionary entries.  Provide that capability
     * to other classes w/i the package.
     * @param contig New value for the contig.
     */
    protected void setContig(String contig) {
        this.contig = contig;
    }

    /**
     * Gets the location of this contig within the fasta.
     * @return seek position within the fasta.
     */
    public long getLocation() {
        return location;
    }

    /**
     * Gets the size, in bytes, of the data in the contig.
     * @return size of the contig bases in bytes.
     */
    public long getSize() {
        return size;
    }

    /**
     * Gets the number of bases in a given line.
     * @return Number of bases in the fasta line.
     */
    public int getBasesPerLine() {
        return basesPerLine;
    }

    /**
     * How many bytes (bases + whitespace) are consumed by the
     * given line?
     * @return Number of bytes in a line.
     */
    public int getBytesPerLine() {
        return bytesPerLine;
    }

    public int getSequenceIndex() {
        return sequenceIndex;
    }

    /**
     * For debugging.  Emit the contents of each contig line.
     * @return A string representation of the contig line.
     */
    public String toString() {
        return String.format("contig %s; location %d; size %d; basesPerLine %d; bytesPerLine %d", contig,
                location,
                size,
                basesPerLine,
                bytesPerLine);
    }

    /**
     * Compare this index entry to another index entry.
     * @param other another FastaSequenceIndexEntry
     * @return True if each has the same name, location, size, basesPerLine and bytesPerLine
     */
    public boolean equals(Object other) {
        if(!(other instanceof HadoopFastaSequenceIndexEntry))
            return false;

        if (this == other) return true;

        HadoopFastaSequenceIndexEntry otherEntry = (HadoopFastaSequenceIndexEntry)other;
        return (contig.equals(otherEntry.contig) && size == otherEntry.size && location == otherEntry.location
                && basesPerLine == otherEntry.basesPerLine && bytesPerLine == otherEntry.bytesPerLine);
    }

    /**
     * In general, we expect one entry per contig, so compute the hash based only on the contig.
     * @return A unique hash code representing this object.
     */
    public int hashCode() {
        return contig.hashCode();
    }
}
