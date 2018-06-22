package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.*;
import java.util.List;

@DefaultSerializer(ContigScorer.Serializer.class)
public class ContigScorer implements Serializable {
    private static final long serialVersionUID = 1L;

    private final int nDivisions;
    private final double stepSplits;
    private final double stepCoverage;
    private final double maxBin; // number of counts in histogram bin having the maximum number of counts as a double
    private final IntHistogram histogram;

    public ContigScorer( final List<AlignedAssemblyOrExcuse> assemblies, final double genomicCoverage ) {
        double maxSplits = 0.;
        double maxCoverage = 0.;
        int nContigs = 0;
        for ( final AlignedAssemblyOrExcuse alignedAssemblyOrExcuse : assemblies ) {
            if ( alignedAssemblyOrExcuse.isNotFailure() ) {
                for ( final ContigScore contigScore : alignedAssemblyOrExcuse.getContigScores() ) {
                    maxSplits = Math.max(maxSplits, contigScore.getSplitReadsPerBase());
                    maxCoverage = Math.max(maxCoverage, contigScore.getMeanCoverage());
                    nContigs += 1;
                }
            }
        }

        // cap coverage -- sometimes there are crazy high values
        if ( maxCoverage > 3.*genomicCoverage ) {
            maxCoverage = 3.*genomicCoverage;
        }

        // an average of 20 counts per division, but no fewer than 2 nor more than 50
        nDivisions = Math.max(50, Math.min(2, nContigs / 20));
        stepSplits = 2.*getStep(maxSplits, nDivisions);
        stepCoverage = getStep(maxCoverage, nDivisions);
        // allocate an extra row for coverage overflow counts
        histogram = new IntHistogram(nDivisions *(nDivisions +1));

        for ( final AlignedAssemblyOrExcuse alignedAssemblyOrExcuse : assemblies ) {
            if ( alignedAssemblyOrExcuse.isNotFailure() ) {
                for ( final ContigScore contigScore : alignedAssemblyOrExcuse.getContigScores() ) {
                    histogram.addObservation(getBin(contigScore));
                }
            }
        }

        maxBin = histogram.getMaxNObservations();
        try ( final BufferedWriter writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile("hdfs://svdev-tws-m:8020/user/tsharpe/wgs1/contigScorer.txt"))) ) {
            writer.write("nContigs = " + nContigs); writer.newLine();
            writer.write( "maxSplits = " + maxSplits); writer.newLine();
            writer.write("maxCoverage = " + maxCoverage); writer.newLine();
            writer.write("nDivisions = " + nDivisions); writer.newLine();
            writer.write("stepSplits = " + stepSplits); writer.newLine();
            writer.write("stepCoverage = " + stepCoverage); writer.newLine();
            writer.write( "maxBin = " + maxBin); writer.newLine();
            for ( int idx = 0; idx != histogram.getMaximumTrackedValue(); ++idx ) {
                writer.write("hist[" + idx + "] = " + histogram.getNObservations(idx)); writer.newLine();
            }
            writer.write( "hist[overflow] = " + histogram.getNObservations(histogram.getMaximumTrackedValue()+1));  writer.newLine();
            writer.write("hist total obs = " + histogram.getTotalObservations()); writer.newLine();
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write ContigScorer state");
        }
    }

    public ContigScorer( final Kryo kryo, final Input input ) {
        nDivisions = input.readInt();
        stepSplits = input.readDouble();
        stepCoverage = input.readDouble();
        maxBin = input.readDouble();
        final IntHistogram.Serializer intHistogramSerializer = new IntHistogram.Serializer();
        histogram = intHistogramSerializer.read(kryo, input, IntHistogram.class);
    }

    public void serialize( final Kryo kryo, final Output output ) {
        output.writeInt(nDivisions);
        output.writeDouble(stepSplits);
        output.writeDouble(stepCoverage);
        output.writeDouble(maxBin);
        final IntHistogram.Serializer intHistogramSerializer = new IntHistogram.Serializer();
        intHistogramSerializer.write(kryo, output, histogram);
    }

    public float score( final ContigScore contigScore ) {
        if ( maxBin == 0. ) {
            throw new GATKException("there are no observations to score against");
        }
        // we don't test for coverage being too high, because getBin has an overflow capability in that dimension
        if ( contigScore.getSplitReadsPerBase() < 0. ||
                contigScore.getSplitReadsPerBase()/ stepSplits >= nDivisions ||
                contigScore.getMeanCoverage() < 0. ) {
            throw new GATKException("contigScore out of bounds");
        }
        return (float)(histogram.getNObservations(getBin(contigScore))/maxBin);
    }

    private static double getStep( final double maxVal, final int nDivisions ) {
        if ( maxVal == 0. ) return 0.5;

        double step = maxVal / nDivisions;
        while ( (int)(maxVal/step) == nDivisions ) {
            step = Math.nextUp(step);
        }
        return step;
    }

    private int getBin( final ContigScore contigScore ) {
        final int column = Math.min((int)(contigScore.getSplitReadsPerBase() / stepSplits), nDivisions-1);
        final int row = Math.min(nDivisions, (int)(contigScore.getMeanCoverage() / stepCoverage));
        return nDivisions * row + column;
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<ContigScorer> {
        @Override
        public void write( final Kryo kryo, final Output output, final ContigScorer contigScorer ) {
            contigScorer.serialize(kryo, output);
        }

        @Override
        public ContigScorer read( final Kryo kryo, final Input input, final Class<ContigScorer> klass ) {
            return new ContigScorer(kryo, input);
        }
    }

    @DefaultSerializer(ContigScore.Serializer.class)
    public static class ContigScore {
        private final double splitReadsPerBase;
        private final double meanCoverage;

        public ContigScore( final double splitReadsPerBase, final double meanCoverage ) {
            this.splitReadsPerBase = splitReadsPerBase;
            this.meanCoverage = meanCoverage;
        }

        public ContigScore( final Kryo kryo, final Input input ) {
            splitReadsPerBase = input.readDouble();
            meanCoverage = input.readDouble();
        }

        public double getSplitReadsPerBase() { return splitReadsPerBase; }
        public double getMeanCoverage() { return meanCoverage; }

        public void serialize( final Kryo kryo, final Output output ) {
            output.writeDouble(splitReadsPerBase);
            output.writeDouble(meanCoverage);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<ContigScore> {
            @Override
            public void write( final Kryo kryo, final Output output, final ContigScore contigScore ) {
                contigScore.serialize(kryo, output);
            }

            @Override
            public ContigScore read( final Kryo kryo, final Input input, final Class<ContigScore> klass ) {
                return new ContigScore(kryo, input);
            }
        }
    }
}
