package htsjdk.samtools.util;


/**
 * Minimal immutable class representing a 1-based closed ended genomic interval
 * SimpleInterval does not allow null contig names.  It cannot represent an unmapped Locatable.
 *
 *@warning 0 length intervals are allowed, which are represented as end = start-1
 */
public class SimpleInterval implements Locatable {

    private final int start;
    private final int end;
    private final String contig;

    /**
     * Create a new immutable 1-based interval of the form [start, end]
     * @param contig the name of the contig, must not be null
     * @param start  1-based inclusive start position
     * @param end  1-based inclusive end position
     */
    public SimpleInterval(final String contig, final int start, final int end){
        if(contig == null){
            throw new IllegalArgumentException("contig cannot be null");
        }
        if (start <= 0 ){
            throw new IllegalArgumentException("SimpleInterval is 1 based, so start must be >= 1, start:%d" + start);
        }
        if (end < start -1){
            throw new IllegalArgumentException( String.format("end must be >= start -1 . start:%d end:%d", start, end));
        }
        this.contig = contig;
        this.start = start;
        this.end = end;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SimpleInterval that = (SimpleInterval) o;

        if (end != that.end) return false;
        if (start != that.start) return false;
        if (!contig.equals(that.contig)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = start;
        result = 31 * result + end;
        result = 31 * result + contig.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return String.format("%s:%s-%s", contig, start, end);
    }

    /**
     * @return name of the contig this is mapped to
     */
    public final String getContig(){
        return contig;
    }

    /** Gets the 1-based start position of the interval on the contig. */
    public final int getStart(){
        return start;
    }

    /**
     * @return the 1-based closed-ended end position of the interval on the contig.
     * @warning in the case of a 0-length interval getEnd() == getStart - 1 */
    public final int getEnd(){
        return end;
    }
}
