package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;

/** statistics of fragment length distributions */
@DefaultSerializer(LibraryStatistics.Serializer.class)
public final class LibraryStatistics {
    private final IntHistogram.CDF fragmentSizeCDF;
    private final int median;
    private final float negativeMAD;
    private final float positiveMAD;
    private int coverage;
    private long nReads;
    private float readStartFrequency; // measured in read starts per reference base

    public LibraryStatistics( final IntHistogram.CDF fragmentSizeCDF,
                              final long nBases, final long nReads, final long nRefBases ) {
        this.fragmentSizeCDF = fragmentSizeCDF;
        median = this.fragmentSizeCDF.median();
        negativeMAD = this.fragmentSizeCDF.leftMedianDeviation(median);
        positiveMAD = this.fragmentSizeCDF.rightMedianDeviation(median);
        coverage = (int)((nBases+nRefBases/2)/nRefBases);
        this.nReads = nReads;
        readStartFrequency = 1.f*nReads/nRefBases;
    }

    private LibraryStatistics( final Kryo kryo, final Input input ) {
        fragmentSizeCDF = new IntHistogram.CDF.Serializer().read(kryo, input, IntHistogram.CDF.class);
        median = fragmentSizeCDF.median();
        negativeMAD = fragmentSizeCDF.leftMedianDeviation(median);
        positiveMAD = fragmentSizeCDF.rightMedianDeviation(median);
        coverage = input.readInt();
        nReads = input.readLong();
        readStartFrequency = input.readFloat();
    }

    private void serialize( final Kryo kryo, final Output output ) {
        new IntHistogram.CDF.Serializer().write(kryo, output, fragmentSizeCDF);
        output.writeInt(coverage);
        output.writeLong(nReads);
        output.writeFloat(readStartFrequency);
    }

    public IntHistogram.CDF getCDF() { return fragmentSizeCDF; }
    public int getMedian() { return median; }
    public float getNegativeMAD() { return negativeMAD; }
    public float getPositiveMAD() { return positiveMAD; }
    public int getCoverage() { return coverage; }
    public long getNReads() { return nReads; }
    public float getReadStartFrequency() { return readStartFrequency; }
    public float getZishScore( final int fragmentSize ) {
        if ( fragmentSize < 0 ) {
            throw new GATKException("negative fragment size");
        }
        final int diff = fragmentSize - median;
        if ( diff == 0 ) return 0.0f;
        if ( diff > 0 ) return 1.0f * diff / positiveMAD;
        return 1.0f * diff / negativeMAD;
    }

    public IntHistogram createEmptyHistogram() { return fragmentSizeCDF.createEmptyHistogram(); }

    public final static class Serializer extends com.esotericsoftware.kryo.Serializer<LibraryStatistics> {
        @Override
        public void write( final Kryo kryo, final Output output, final LibraryStatistics stats ) {
            stats.serialize(kryo, output);
        }

        @Override
        public LibraryStatistics read( final Kryo kryo, final Input input,
                                       final Class<LibraryStatistics> klass ) {
            return new LibraryStatistics(kryo, input);
        }
    }
}
