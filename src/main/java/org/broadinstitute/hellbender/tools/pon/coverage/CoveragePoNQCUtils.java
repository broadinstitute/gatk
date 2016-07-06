package org.broadinstitute.hellbender.tools.pon.coverage;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

/**
 * Class for doing QC on a PoN
 */
public final class CoveragePoNQCUtils {
    private static final Logger logger = LogManager.getLogger(CoveragePoNQCUtils.class);

    static final double AMP_THRESHOLD = 1.12;
    static final double DEL_THRESHOLD = 0.88;

    private CoveragePoNQCUtils() {}

    /**
     * Return list of samples in the given PoN that likely contain arm level events.  Samples detected will have to be
     *  removed.
     *
     * @param qcPoN -- a QC PoN, which typically has fewer eigensamples.
     * @param ctx -- use null if not using Spark
     * @return never {@code null}
     */
    public static <T extends CoveragePoNNormalizationResult> List<String> retrieveSamplesWithArmLevelEvents(final CoveragePanelOfNormals<T> qcPoN, final JavaSparkContext ctx) {
        Utils.nonNull(qcPoN, "QC PoN cannot be null.");
        // For each sample in the PoN, tangent normalize against the qc reduced PoN.
        logger.info("Tangent normalizing all normals in the PoN...");
        final CoveragePoNNormalizationResult qcNormalization = qcPoN.normalizeNormalsInPoN(ctx);

        // Get the median target CR estimate per contig across all samples.
        final Map<String, Double> contigToMedian = getContigToMedianCRMap(qcNormalization.getProfile());

        // Segment the result for each tangent normalization result from the QC PoN
        if (ctx == null) {
            final List<ReadCountCollection> singleSampleTangentNormalizedReadCounts = createIndividualReadCountCollections(qcNormalization.getProfile());
            return identifySamplesWithSuspiciousContigs(singleSampleTangentNormalizedReadCounts, contigToMedian);
        }
        final JavaRDD<ReadCountCollection> parallelSingleSampleTangentNormalizedReadCounts = createParallelIndividualReadCountCollections(qcNormalization.getProfile(), ctx);
        return identifySamplesWithSuspiciousContigs(parallelSingleSampleTangentNormalizedReadCounts, ctx, contigToMedian);
    }

    @VisibleForTesting
    static Map<String, Double> getContigToMedianCRMap(final ReadCountCollection readCountCollection) {
        final List<String> allContigsPresent = retrieveAllContigsPresent(readCountCollection);
        final Map<String, Double> contigToMedian = new LinkedHashMap<>();
        for (String contig: allContigsPresent) {
            final ReadCountCollection oneContigReadCountCollection = readCountCollection.subsetTargets(readCountCollection.targets().stream().filter(t -> t.getContig().equals(contig)).collect(Collectors.toSet()));
            final double[] flatCounts = Doubles.concat(oneContigReadCountCollection.counts().getData());

            // Put into CRSpace
            final double[] flatCountsInCRSpace = DoubleStream.of(flatCounts).map(d -> Math.pow(2, d)).toArray();
            contigToMedian.put(contig, new Median().evaluate(flatCountsInCRSpace));
        }
        return contigToMedian;
    }

    /**
     *  Given a single sample tangent normalization (or other coverage profile), determine whether any contig looks like
     *   it has an arm level event (defined as 25% (or more) of the contig amplified/deleted)
     *
     * @param singleSampleTangentNormalized Tangent normalized data for a single sample.
     * @return never {@code null}
     */
    private static Boolean hasSuspiciousContigs(final ReadCountCollection singleSampleTangentNormalized, final Map<String, Double> contigToMedian) {
        final List<String> allContigsPresent = retrieveAllContigsPresent(singleSampleTangentNormalized);
        for (String contig: allContigsPresent) {
            final ReadCountCollection oneContigReadCountCollection = singleSampleTangentNormalized.subsetTargets(singleSampleTangentNormalized.targets().stream().filter(t -> t.getContig().equals(contig)).collect(Collectors.toSet()));
            final RealVector counts = oneContigReadCountCollection.counts().getColumnVector(0);
            for (int i = 0; i < 4; i++) {
                final RealVector partitionCounts = counts.getSubVector(i * counts.getDimension() / 4, counts.getDimension() / 4);
                final double[] partitionArray = DoubleStream.of(partitionCounts.toArray()).map(d -> Math.pow(2, d)).sorted().toArray();
                double median = new Median().evaluate(partitionArray);
                final double medianShiftInCRSpace = contigToMedian.getOrDefault(contig, 1.0) - 1.0;
                median -= medianShiftInCRSpace;
                if ((median > AMP_THRESHOLD) || (median < DEL_THRESHOLD)) {
                    logger.info("Suspicious contig: " + singleSampleTangentNormalized.columnNames().get(0) + " " + contig + " (" + median + " -- " + i + ")");
                    return true;
                }
            }
        }
        return false;
    }

    @VisibleForTesting
    static List<String> identifySamplesWithSuspiciousContigs(final List<ReadCountCollection> singleSampleTangentNormalizedReadCounts, final Map<String, Double> contigToMedian) {
        logger.info("Checking for samples with suspicious contigs...");
        return singleSampleTangentNormalizedReadCounts.stream().filter(rcc -> hasSuspiciousContigs(rcc, contigToMedian)).map(rcc -> rcc.columnNames().get(0)).collect(Collectors.toList());
    }

    @VisibleForTesting
    static List<String> identifySamplesWithSuspiciousContigs(final JavaRDD<ReadCountCollection> parallelReadCountCollections, final JavaSparkContext ctx, final Map<String, Double> contigToMedian) {
        Broadcast<Map<String, Double>> broadcastedContigToMedian = ctx.broadcast(contigToMedian);
        return parallelReadCountCollections.filter(rcc -> hasSuspiciousContigs(rcc, broadcastedContigToMedian.value())).map(rcc -> rcc.columnNames().get(0)).collect();
    }

    /**
     *  Split a read count collection into a separate read count collection for each sample
     *
     * @param readCountCollection -- input read count collection assumed to have 1+ samples
     * @return not {@code null}
     */
    @VisibleForTesting
    static List<ReadCountCollection> createIndividualReadCountCollections(final ReadCountCollection readCountCollection) {
        logger.info("Convert ReadCountCollection with multiple samples into list of read count collections each with one sample...");
        final List<String> sampleNames = readCountCollection.columnNames();
        return sampleNames.stream().map(s -> readCountCollection.subsetColumns(Collections.singleton(s))).collect(Collectors.toList());
    }

    /**
     *  Split a read count collection into a separate read count collection for each sample.
     *
     * @param readCountCollection -- input read count collection assumed to have 1+ samples
     * @param ctx Use {@code null} if no spark context is available.  Result serialization is the big bottleneck, so it may be worthwhile to skip spark with fewer available cores
     * @return not {@code null}
     */
    @VisibleForTesting
    static JavaRDD<ReadCountCollection> createParallelIndividualReadCountCollections(final ReadCountCollection readCountCollection, final JavaSparkContext ctx) {
        final List<String> sampleNames = readCountCollection.columnNames();
        Broadcast<ReadCountCollection> broadcastedReadCountCollection = ctx.broadcast(readCountCollection);
        JavaRDD<String> parallelSampleNames = ctx.parallelize(sampleNames, Math.max(sampleNames.size()/10, 4));
        return parallelSampleNames.map(s -> broadcastedReadCountCollection.value().subsetColumns(Collections.singleton(s)));
    }

    private static List<String> retrieveAllContigsPresent(final ReadCountCollection readCountCollection) {
        return readCountCollection.targets().stream().map(Target::getContig).distinct().collect(Collectors.toList());
    }
}
