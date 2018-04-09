package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.FermiLiteAssemblyHandler.ContigScore;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;

import java.io.Serializable;
import java.util.List;

@DefaultSerializer(ContigScorer.Serializer.class)
public class ContigScorer implements Serializable {
    private static final long serialVersionUID = 1L;

    private final IntHistogram meanASHistogram;
    private final IntHistogram meanCoverageHistogram;
    private final float maxASBin;
    private final float maxCoverageBin;

    public ContigScorer( final List<AlignedAssemblyOrExcuse> assemblies ) {
        float maxMeanAS = 0.f;
        float maxMeanCoverage = 0.f;
        for ( final AlignedAssemblyOrExcuse alignedAssemblyOrExcuse : assemblies ) {
            if ( alignedAssemblyOrExcuse.isNotFailure() ) {
                for ( final ContigScore contigScore : alignedAssemblyOrExcuse.getContigScores() ) {
                    if ( contigScore.getMeanAS() > maxMeanAS ) {
                        maxMeanAS = contigScore.getMeanAS();
                    }
                    if ( contigScore.getMeanCoverage() > maxMeanCoverage ) {
                        maxMeanCoverage = contigScore.getMeanCoverage();
                    }
                }
            }
        }

        meanASHistogram = new IntHistogram((int)maxMeanAS);
        meanCoverageHistogram = new IntHistogram((int)maxMeanCoverage);

        for ( final AlignedAssemblyOrExcuse alignedAssemblyOrExcuse : assemblies ) {
            if ( alignedAssemblyOrExcuse.isNotFailure() ) {
                for ( final ContigScore contigScore : alignedAssemblyOrExcuse.getContigScores() ) {
                    meanASHistogram.addObservation((int)contigScore.getMeanAS());
                    meanCoverageHistogram.addObservation((int)contigScore.getMeanCoverage());
                }
            }
        }

        maxASBin = (float)meanASHistogram.getMaxNObservations();
        maxCoverageBin = (float)meanCoverageHistogram.getMaxNObservations();
    }

    public ContigScorer( final Kryo kryo, final Input input ) {
        final IntHistogram.Serializer intHistogramSerializer = new IntHistogram.Serializer();
        meanASHistogram = intHistogramSerializer.read(kryo, input, IntHistogram.class);
        meanCoverageHistogram = intHistogramSerializer.read(kryo, input, IntHistogram.class);
        maxASBin = input.readFloat();
        maxCoverageBin = input.readFloat();
    }

    public void serialize( final Kryo kryo, final Output output ) {
        final IntHistogram.Serializer intHistogramSerializer = new IntHistogram.Serializer();
        intHistogramSerializer.write(kryo, output, meanASHistogram);
        intHistogramSerializer.write(kryo, output, meanCoverageHistogram);
        output.writeFloat(maxASBin);
        output.writeFloat(maxCoverageBin);
    }

    public float score( final ContigScore score ) {
        final float asVal = meanASHistogram.getNObservations((int)score.getMeanAS()) / maxASBin;
        final float coverageVal = meanCoverageHistogram.getNObservations((int)score.getMeanCoverage()) / maxCoverageBin;
        return asVal*coverageVal;
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
}
