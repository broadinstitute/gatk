package org.broadinstitute.hellbender.fakedata;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.tools.exome.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * GC bias data consists of a {@link ReadCountCollection} and an array of gc contents
 * in order of the read counts' targets.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class GCBiasSimulatedData {

    // a reasonable default GC bias curve
    private static Function<Double, Double> QUADRATIC_GC_BIAS_CURVE = gc -> 0.5 + 2*gc*(1-gc);

    // visible for the integration test
    public static Pair<ReadCountCollection, double[]> simulatedData(final int numTargets, final int numSamples) {
        final List<Target> phonyTargets = SimulatedTargets.phonyTargets(numTargets);
        final List<String> phonySamples = SimulatedSamples.phonySamples(numSamples);

        final Random random = new Random(13);
        final double[] gcContentByTarget = IntStream.range(0, numTargets)
                .mapToDouble(n -> 0.5 + 0.2*random.nextGaussian())
                .map(x -> Math.min(x,0.95)).map(x -> Math.max(x,0.05)).toArray();
        final double[] gcBiasByTarget = Arrays.stream(gcContentByTarget).map(QUADRATIC_GC_BIAS_CURVE::apply).toArray();

        // model mainly GC bias with a small random amount of non-GC bias
        // thus noise after GC correction should be nearly zero
        final RealMatrix counts = new Array2DRowRealMatrix(numTargets, numSamples);
        counts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int target, final int column, final double value) {
                return gcBiasByTarget[target]*(1.0 + 0.01*random.nextDouble());
            }
        });
        final ReadCountCollection rcc = new ReadCountCollection(phonyTargets, phonySamples, counts);
        return new ImmutablePair<>(rcc, gcContentByTarget);
    }

    /**
     *
     * @param readCountsFile    A simulated read counts file with GC bias effects
     * @param targetsFile       A simulated targets file with GC content annotation
     */
    public static void makeGCBiasInputFiles(final Pair<ReadCountCollection, double[]> data,
                                            final File readCountsFile, final File targetsFile) throws IOException {
        final ReadCountCollection inputCounts = data.getLeft();
        final double[] gcContentByTarget = data.getRight();
        ReadCountCollectionUtils.write(readCountsFile, inputCounts);

        final TargetWriter writer = new TargetWriter(targetsFile, Collections.singleton(TargetAnnotation.GC_CONTENT));
        for (int i = 0; i < gcContentByTarget.length; i++) {
            final Target unannotatedTarget = inputCounts.records().get(i).getTarget();
            final TargetAnnotationCollection annotations = new TargetAnnotationCollection();
            annotations.put(TargetAnnotation.GC_CONTENT, Double.toString(gcContentByTarget[i]));
            final Target target = new Target(unannotatedTarget.getName(), unannotatedTarget.getInterval(), annotations);
            writer.writeRecord(target);
        }
        writer.close();
    }
}
