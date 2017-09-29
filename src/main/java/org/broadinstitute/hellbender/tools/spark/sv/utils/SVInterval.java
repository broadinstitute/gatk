package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Naturally collating, simple interval.
 * Some methods assume that the interval is half-open.
 */
@DefaultSerializer(SVInterval.Serializer.class)
@VisibleForTesting
public final class SVInterval implements Comparable<SVInterval> {
    private final int contig;
    private final int start;
    private final int end;

    public SVInterval( final int contig, final int start, final int end ) {
        this.contig = contig;
        this.start = start;
        this.end = end;
    }

    private SVInterval( final Kryo kryo, final Input input ) {
        contig = input.readInt();
        start = input.readInt();
        end = input.readInt();
    }

    private void serialize( final Kryo kryo, final Output output ) {
        output.writeInt(contig);
        output.writeInt(start);
        output.writeInt(end);
    }

    public int getContig() { return contig; }
    public int getStart() { return start; }
    public int getEnd() { return end; }
    public int getLength() { return end-start; }

    /** This definition is appropriate for half-open intervals.
     *  If you're building your intervals as closed, you're going to have trouble.
     */
    public boolean overlaps( final SVInterval that ) {
        return this.contig == that.contig && this.start < that.end && that.start < this.end;
    }

    public boolean isDisjointFrom( final SVInterval that ) {
        return !overlaps(that);
    }

    public boolean isUpstreamOf( final SVInterval that ) {
        return this.contig < that.contig || (this.contig == that.contig && this.end <= that.start);
    }

    public int gapLen( final SVInterval that ) {
        if ( this.contig != that.contig ) return Integer.MAX_VALUE;
        return that.start - this.end;
    }

    public int overlapLen( final SVInterval that ) {
        if ( this.contig != that.contig ) return 0;
        return Math.max(0, Math.min(this.end, that.end) - Math.max(this.start, that.start));
    }

    public SVInterval join( final SVInterval that ) {
        if ( this.contig != that.contig ) throw new GATKException("Joining across contigs.");
        return new SVInterval(contig, Math.min(this.start, that.start), Math.max(this.end, that.end));
    }

    @Override
    public boolean equals( final Object obj ) {
        if ( !(obj instanceof SVInterval) ) return false;
        final SVInterval that = (SVInterval)obj;
        return this.contig == that.contig && this.start == that.start && this.end == that.end;
    }

    @Override
    public int hashCode() {
        return 47*(47*(47*(47*contig)+start)+end);
    }

    @Override
    public int compareTo( final SVInterval that ) {
        int result = Integer.compare(this.contig, that.contig);
        if ( result == 0 ) {
            result = Integer.compare(this.start, that.start);
            if ( result == 0 ) result = Integer.compare(this.end, that.end);
        }
        return result;
    }

    @Override
    public String toString() { return Integer.toString(contig)+"["+start+":"+end+"]"; }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SVInterval> {
        @Override
        public void write( final Kryo kryo, final Output output, final SVInterval interval ) {
            interval.serialize(kryo, output);
        }

        @Override
        public SVInterval read(final Kryo kryo, final Input input, final Class<SVInterval> klass ) {
            return new SVInterval(kryo, input);
        }
    }

}
