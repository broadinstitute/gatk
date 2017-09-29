package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;

/** statistics of fragment length distributions */
@DefaultSerializer(FragmentLengthStatistics.Serializer.class)
public final class FragmentLengthStatistics {
    private final IntHistogram.CDF cdf;
    private final int median;
    private final float negativeMAD;
    private final float positiveMAD;

    public FragmentLengthStatistics( final IntHistogram allLibraryReadsHistogram ) {
        cdf = allLibraryReadsHistogram.getCDF();
        median = cdf.median();
        negativeMAD = cdf.leftMedianDeviation(median);
        positiveMAD = cdf.rightMedianDeviation(median);
    }

    private FragmentLengthStatistics( final Kryo kryo, final Input input ) {
        cdf = new IntHistogram.CDF.Serializer().read(kryo, input, IntHistogram.CDF.class);
        median = cdf.median();
        negativeMAD = cdf.leftMedianDeviation(median);
        positiveMAD = cdf.rightMedianDeviation(median);
    }

    private void serialize( final Kryo kryo, final Output output ) {
        new IntHistogram.CDF.Serializer().write(kryo, output, cdf);
    }

    public IntHistogram.CDF getCDF() { return cdf; }
    public int getMedian() { return median; }
    public float getNegativeMAD() { return negativeMAD; }
    public float getPositiveMAD() { return positiveMAD; }

    public float getZishScore( final int fragmentSize ) {
        if ( fragmentSize < 0 ) {
            throw new GATKException("negative fragment size");
        }
        final int diff = fragmentSize - median;
        if ( diff == 0 ) return 0.0f;
        if ( diff > 0 ) return 1.0f * diff / positiveMAD;
        return 1.0f * diff / negativeMAD;
    }

    IntHistogram createEmptyHistogram() { return cdf.createEmptyHistogram(); }

    boolean isDifferentByKSStatistic( final IntHistogram sampleHistogram, final float significance ) {
        return cdf.isDifferentByKSStatistic(sampleHistogram, significance);
    }

    public final static class Serializer extends com.esotericsoftware.kryo.Serializer<FragmentLengthStatistics> {
        @Override
        public void write( final Kryo kryo, final Output output, final FragmentLengthStatistics stats ) {
            stats.serialize(kryo, output);
        }

        @Override
        public FragmentLengthStatistics read( final Kryo kryo, final Input input,
                                              final Class<FragmentLengthStatistics> klass ) {
            return new FragmentLengthStatistics(kryo, input);
        }
    }
}
