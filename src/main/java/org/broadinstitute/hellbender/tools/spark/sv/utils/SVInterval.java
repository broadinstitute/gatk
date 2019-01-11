package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Naturally collating, simple interval
 *
 * WARNING: THIS IS NOT THE SAME AS THE BED COORDINATE SYSTEM OR {@link SimpleInterval} !!!!!
 *
 * <p>
 *     SVIntervals are comparable, and provide the same total order as coordinate-sorted BAMs.
 *     They can encode any contiguous stretch of bases on the reference,
 *     or strand-sensitively on its reverse complement.
 *     They can encode 0-length intervals (i.e., locations) to indicate things like the precise,
 *     unambiguous location of an insert between two reference bases.
 *     They're super light-weight, and cache friendly (no references),
 *     and they can be compared and tested for equality quickly and locally.
 *     They know how to serialize themselves with Kryo.
 * </p>
 *
 * In general, flexibility is provided by SVInterval, and
 * one should be responsible for carefully being self-consistent.
 *
 * So please refrain from putting methods here that
 * makes any assumptions about coordinate system, 1-based or 0-based.
 *
 * <p>
 *     Some methods assume that the interval is half-open.
 *     Please read the documentations clearly before picking a convention and starting to construct the intervals.
 *     And one should stick to this convention.
 * </p>
 */
@DefaultSerializer(SVInterval.Serializer.class)
@VisibleForTesting
public final class SVInterval implements Comparable<SVInterval> {
    private final int contig;
    private final int start;
    private final int end;

    /**
     * This constructor uses the {@link SVIntervalConstructorArgsValidator#RefuseNegativeContigAndCoordinates} as the validator.
     */
    public SVInterval( final int contig, final int start, final int end ) {
        this(contig, start, end, SVIntervalConstructorArgsValidator.RefuseNegativeContigAndCoordinates);
    }

    @FunctionalInterface
    public interface SVIntervalConstructorArgsValidator {
        void accept(final int contig, final int start, final int end);

        default SVIntervalConstructorArgsValidator andThen(final SVIntervalConstructorArgsValidator after) {
            Utils.nonNull(after);
            return (final int contig, final int start, final int end) -> {
                accept(contig, start, end); after.accept(contig, start, end);
            };
        }

        SVIntervalConstructorArgsValidator ACCEPTS_ALL = ((contig, start, end) -> {});

        SVIntervalConstructorArgsValidator RefuseNegativeContigs =
                ((contig, start, end) -> {
                    if (contig < 0)
                        throw new IllegalArgumentException("provided contig is negative: " + contig);
                });

        SVIntervalConstructorArgsValidator RefuseNegativeCoordinates =
                ((contig, start, end) -> {
                    if (start < 0)
                        throw new IllegalArgumentException("provided start is negative: " + start);
                    if (end < 0)
                        throw new IllegalArgumentException("provided end is negative: " + end);
                });

        SVIntervalConstructorArgsValidator RefuseNegativeContigAndCoordinates =
                RefuseNegativeContigs.andThen(RefuseNegativeCoordinates);
    }

    public SVInterval( final int contig, final int start, final int end, final SVIntervalConstructorArgsValidator argsValidator) {
        argsValidator.accept(contig, start, end);
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

    // assumes the interval is half-open
    public int getLength() { return end - start; }

    public SVLocation getStartLocation() { return new SVLocation(contig, start); }

    /**
     * This definition is appropriate for half-open intervals.
     * If you're building your intervals as closed, you're going to have trouble.
     */
    public boolean overlaps( final SVInterval that ) {
        return this.contig == that.contig && this.start < that.end && that.start < this.end;
    }

    // assumes the interval is half-open
    public boolean isDisjointFrom( final SVInterval that ) {
        return !overlaps(that);
    }

    /**
     * Returns false if the two intervals overlap as judged by {@link #overlaps(SVInterval)}
     */
    public boolean isUpstreamOf( final SVInterval that ) {
        return this.contig < that.contig || (this.contig == that.contig && this.end <= that.start);
    }

    public boolean contains( final SVInterval that ) {
        return this.contig == that.contig && this.start <= that.start && this.end >= that.end;
    }

    /**
     * Assumes interval half-open and this {@link #isUpstreamOf(SVInterval)}.
     */
    public int gapLen( final SVInterval that ) {
        if ( this.contig != that.contig ) return Integer.MAX_VALUE;
        return that.start - this.end;
    }

    // assumes half-open
    public int overlapLen( final SVInterval that ) {
        if ( this.contig != that.contig ) return 0;
        return Math.max(0, Math.min(this.end, that.end) - Math.max(this.start, that.start));
    }

    // no assumption about half-open
    public SVInterval join( final SVInterval that ) {
        if ( this.contig != that.contig ) throw new GATKException("Joining across contigs.");
        return new SVInterval(contig, Math.min(this.start, that.start), Math.max(this.end, that.end));
    }

    /**
     * Returns {@code null} if the two intervals don't overlap as judged by {@link #overlaps(SVInterval)}.
     */
    public SVInterval intersect(final SVInterval that) {
        if (! this.overlaps(that)) return null;
        return new SVInterval(this.getContig(), Math.max(this.start, that.start), Math.min(this.end, that.end));
    }

    // a floored mid point between start and end, e.g. if start == 1 and end == 2, the midpoint will be (1+2)/2 -> 1.
    public int midpoint() {
        return (start + end) / 2;
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
