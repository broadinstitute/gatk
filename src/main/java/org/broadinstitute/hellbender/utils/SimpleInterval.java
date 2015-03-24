package org.broadinstitute.hellbender.utils;


import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Minimal immutable class representing a 1-based closed ended genomic interval
 * SimpleInterval does not allow null contig names.  It cannot represent an unmapped Locatable.
 *
 *@warning 0 length intervals are NOT currently allowed, but support may be added in the future
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
        if ( contig == null ) {
            throw new IllegalArgumentException("contig cannot be null");
        }
        if ( start <= 0 ) {
            throw new IllegalArgumentException("SimpleInterval is 1 based, so start must be >= 1, start:%d" + start);
        }
        if ( end < start ) {
            throw new IllegalArgumentException(String.format("end must be >= start. start:%d end:%d", start, end));
        }
        this.contig = contig;
        this.start = start;
        this.end = end;
    }

    /**
     * HACK to avoid refactoring the interval parsing system for now.
     * Create a SimpleInterval from a String that has the format chr:start-end
     * @param intervalString in the format chr:start-end
     */
    public static SimpleInterval valueOf(String intervalString){
        String[] pieces = intervalString.split("[:-]");
        try{
            if( pieces.length != 3){
                throw new GATKException("Failed to create a new SimpleInterval from the value, expected a string of the form <chr>:<start>-<end>, got %s.");
            }
            return new SimpleInterval(pieces[0], Integer.valueOf(pieces[1]), Integer.valueOf(pieces[2]));
        } catch( NumberFormatException e ) {
            throw new GATKException(String.format("Failed to create a new SimpleInterval from the value, expected a string of the form <chr>:<start>-<end>, got %s.",
                    intervalString), e);
        }
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
     */
    public final int getEnd(){
        return end;
    }

    /**
     * @return number of bases covered by this interval (will always be > 0)
     */
    public final int size() {
        return end - start + 1;
    }

    /**
     * Determines whether this interval overlaps the provided interval.
     *
     * @param other interval to check
     * @return true if this interval overlaps other, otherwise false
     */
    public final boolean overlaps( final SimpleInterval other ) {
        if ( other == null ) {
            return false;
        }

        return this.contig.equals(other.contig) && this.start <= other.end && other.start <= this.end;
    }

    /**
     * Determines whether this interval contains the entire region represented by other
     *
     * @param other interval to check
     * @return true if this interval contains all of the bases spanned by other, otherwise false
     */
    public final boolean contains( final SimpleInterval other ) {
        if ( other == null ) {
            return false;
        }

        return this.contig.equals(other.contig) && this.start <= other.start && this.end >= other.end;
    }
}
