package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.Utils;
import picard.sam.util.Pair;

import java.util.*;

public class DepthMatrixLoader {

    protected static final int MIN_BINS_FOR_TRIMMING = 4;

    protected final FeatureDataSource<DepthEvidence> dataSource;
    protected final SAMSequenceDictionary dictionary;
    protected final int bins;
    protected final long largeSizeCutoff;
    protected final int largeVariantRegionPoints;
    protected final int largeVariantWindow;

    protected record CompressionResult(double compression, SimpleInterval adjustedRegion) {
    }

    public DepthMatrixLoader(final FeatureDataSource<DepthEvidence> dataSource,
                             final int bins, final long largeSizeCutoff, final int largeVariantRegionPoints, final int largeVariantWindow,
                             final SAMSequenceDictionary dictionary) {
        Utils.nonNull(dataSource);
        Utils.validateArg(bins > 0, "bins must be greater than zero");
        Utils.validateArg(largeSizeCutoff >= 0, "large size cutoff must be non-negative");
        Utils.validateArg(largeVariantRegionPoints > 0, "large variant region points must be greater than zero");
        Utils.validateArg(largeVariantWindow >= 0, "large variant region window must be non-negative");
        Utils.nonNull(dictionary);
        this.dataSource = dataSource;
        this.dictionary = dictionary;
        this.bins = bins;
        this.largeSizeCutoff = largeSizeCutoff;
        this.largeVariantRegionPoints = largeVariantRegionPoints;
        this.largeVariantWindow = largeVariantWindow;
    }

    public DepthMatrix load(final SimpleInterval interval, final Map<String, Double> sampleMedians) {
        Utils.nonNull(interval);
        if (!IntervalUtils.intervalIsOnDictionaryContig(interval, dictionary)) {
            throw new GATKException("Load interval " + interval + " failed to validate against the sequence dictionary");
        }
        Utils.nonNull(sampleMedians);

        final DepthMatrix rawMatrix;
        if (interval.getLengthOnReference() > largeSizeCutoff) {
            rawMatrix = DepthMatrix.fromDataSourceSubsampled(interval, dataSource, largeVariantWindow, largeVariantRegionPoints, dictionary);
        } else {
            rawMatrix = DepthMatrix.fromDataSource(interval, dataSource);
        }

        // Calculate compression and region bounds
        final CompressionResult compressionResult = calculateCompression(rawMatrix, interval, bins);

        // Trim bins from compression
        final DepthMatrix trimmedMatrix = trimMatrixBeforeCompression(rawMatrix, compressionResult);

        // Set 0 entries to 1
        setZeroesToOnes(trimmedMatrix);

        // Apply compression and normalization
        final DepthMatrix compressedMatrix = compressMatrix(trimmedMatrix, compressionResult);
        final Map<String, Double> compressedMedians = compressMedians(sampleMedians, compressionResult);

        // Normalize by sample depth
        final DepthMatrix normalizedMatrix = normalizeMatrix(compressedMatrix, compressedMedians);
        return trimOuterBins(normalizedMatrix);
    }

    protected static void setZeroesToOnes(final DepthMatrix matrix) {
        // Replace zeros with ones
        for (String sample : matrix.getSampleSet()) {
            final double[] values = matrix.getSample(sample);
            for (int i = 0; i < values.length; i++) {
                if (values[i] == 0.0) {
                    values[i] = 1.0;
                }
            }
        }
    }

    protected static DepthMatrix trimMatrixBeforeCompression(final DepthMatrix matrix, final CompressionResult compressionResult) {
        if (compressionResult == null) {
            return matrix;
        }
        final int trimStart = compressionResult.adjustedRegion.getStart();
        final int trimEnd = compressionResult.adjustedRegion.getEnd();
        final Set<Integer> indicesToRemove = new HashSet<>();
        final List<SimpleInterval> matrixBins = matrix.getBins();
        for (int i = 0; i < matrixBins.size(); i++) {
            final SimpleInterval bin = matrixBins.get(i);
            if (bin.getEnd() <= trimStart || bin.getStart() >= trimEnd) {
                indicesToRemove.add(i);
            }
        }
        if (indicesToRemove.isEmpty()) {
            return matrix;
        }

        final int newNumBins = matrixBins.size() - indicesToRemove.size();
        final List<SimpleInterval> newBins = new ArrayList<>(newNumBins);
        for (int i = 0; i < matrixBins.size(); i++) {
            if (!indicesToRemove.contains(i)) {
                newBins.add(matrixBins.get(i));
            }
        }
        final Map<String, double[]> newMatrix = new HashMap<>();
        for (final String sample : matrix.getSampleSet()) {
            final double[] values = matrix.getSample(sample);
            final double[] newValues = new double[newNumBins];
            int j = 0;
            for (int i = 0; i < values.length; i++) {
                if (!indicesToRemove.contains(i)) {
                    newValues[j++] = values[i];
                }
            }
            newMatrix.put(sample, newValues);
        }
        return new DepthMatrix(newBins, newMatrix);
    }

    protected static DepthMatrix trimOuterBins(final DepthMatrix matrix) {
        final int numBins = matrix.getNumBins();
        if (numBins < MIN_BINS_FOR_TRIMMING) {
            return matrix;
        }
        final List<SimpleInterval> newBins = matrix.getBins().subList(1, numBins - 1);
        final Map<String, double[]> newMatrix = new HashMap<>();
        for (final String sample : matrix.getSampleSet()) {
            final double[] counts = matrix.getSample(sample);
            final double[] trimmedCounts = Arrays.copyOfRange(counts, 1, numBins - 1);
            newMatrix.put(sample, trimmedCounts);
        }
        return new DepthMatrix(newBins, newMatrix);
    }

    protected static CompressionResult calculateCompression(final DepthMatrix matrix, final SimpleInterval interval, final int bins) {
        final int start = interval.getStart() - 1; // covert to 0-based to match R code
        final int end = interval.getEnd();

        // Find bins that overlap with the query interval
        int startBinIdx = -1;
        int endBinIdx = -1;
        for (int i = 0; i < matrix.getBins().size(); i++) {
            final SimpleInterval bin = matrix.getBins().get(i);
            if (startBinIdx == -1 && bin.getStart() - 1 >= start) {
                startBinIdx = i;
            }
            if (bin.getEnd() <= end) {
                endBinIdx = i;
            }
        }

        if (startBinIdx == -1) startBinIdx = 0;
        if (endBinIdx == -1) endBinIdx = matrix.getBins().size() - 1;
        final int numInternalBins = endBinIdx - startBinIdx + 1;

        final int numBins;
        if (numInternalBins < bins) {
            numBins = numInternalBins;
            if (numBins == 0) {
                return null;
            } else if (numBins == 1) {
                final SimpleInterval adjustedRegion = new SimpleInterval(
                        interval.getContig(),
                        matrix.getBins().get(0).getStart(),
                        matrix.getBins().get(matrix.getBins().size() - 1).getEnd()
                );
                return new CompressionResult(1, adjustedRegion);
            }
        } else {
            numBins = bins;
        }

        // Adjust for even compression
        final int remainderForRemoval = (int) Math.floor(((double) (numInternalBins % numBins)) / 2.0);
        final int newStartBinIdx = startBinIdx + remainderForRemoval;
        final int newEndBinIdx = endBinIdx - remainderForRemoval;

        // Calculate compression factor
        final double compression = (newEndBinIdx - newStartBinIdx + 1) / (double) numBins;

        final SimpleInterval adjustedRegion = new SimpleInterval(
                interval.getContig(),
                matrix.getBins().get(newStartBinIdx).getStart(),
                matrix.getBins().get(newEndBinIdx).getEnd()
        );

        return new CompressionResult(compression, adjustedRegion);
    }

    protected static DepthMatrix compressMatrix(final DepthMatrix matrix,
                                                final CompressionResult compressionResult) {
        if (compressionResult == null || compressionResult.compression() <= 1 || matrix.isEmpty()) {
            return matrix;
        }
        final String firstSample = matrix.getSampleSet().iterator().next();
        final List<Pair<Integer, Integer>> compressionIndices = computeCompressionIndices(matrix.getSample(firstSample), compressionResult.compression);
        final List<SimpleInterval> compressedBins = compressBins(matrix.getBins(), compressionIndices);
        final Map<String, double[]> compressedMatrix = new HashMap<>();
        for (final String sample : matrix.getSampleSet()) {
            final double[] counts = matrix.getSample(sample);
            // Apply compression function to this sample's counts
            final double[] compressed = compressCounts(counts, compressionIndices);
            compressedMatrix.put(sample, compressed);
        }
        return new DepthMatrix(compressedBins, compressedMatrix);
    }

    private static Map<String, Double> compressMedians(final Map<String, Double> sampleMedians,
                                                       final CompressionResult compressionResult) {
        if (compressionResult == null || compressionResult.compression() <= 1) {
            return sampleMedians;
        }
        // Approximate rebinned per-sample medians
        final Map<String, Double> adjustedMedians = new HashMap<>();
        for (final String sample : sampleMedians.keySet()) {
            double median = sampleMedians.get(sample);
            if (median == 0.0) median = 1.0;
            adjustedMedians.put(sample, median * compressionResult.compression);
        }
        return adjustedMedians;
    }


    /**
     * Compression function equivalent to the R function
     *
     * @param vals        array of values to compress
     * @param compression compression factor (as double)
     * @return compressed array as double[]
     */
    private static List<Pair<Integer, Integer>> computeCompressionIndices(final double[] vals, final double compression) {
        Utils.nonNull(vals);
        Utils.validateArg(compression > 1., "Compression must be greater than 1");
        if (vals.length == 0) return Collections.emptyList();
        final List<Pair<Integer, Integer>> result = new ArrayList<>();

        // Convert compression to int for array indexing
        final int firstEnd = Math.min(((int) compression) - 1, vals.length - 1);
        result.add(new Pair<>(0, firstEnd));

        // Process remaining chunks
        final int numChunks = (int) (vals.length / compression);
        for (int i = 2; i <= numChunks; i++) {
            final double startIndexFloat = ((i - 1) * compression);
            final int startIndex = (int) startIndexFloat;
            final double endIndexFloat = (i * compression) - 1;
            int endIndex = (int) endIndexFloat;
            // This is a weird edge case to match R indexing behavior
            if (startIndexFloat > startIndex && endIndexFloat == endIndex) {
                endIndex--;
            }
            result.add(new Pair<>(startIndex, endIndex));
        }
        return Collections.unmodifiableList(result);
    }

    private static List<SimpleInterval> compressBins(final List<SimpleInterval> bins, final List<Pair<Integer, Integer>> indices) {
        final List<SimpleInterval> compressedBins = new ArrayList<>(indices.size());
        for (final Pair<Integer, Integer> pair : indices) {
            compressedBins.add(makeCompressedBin(bins, pair.getLeft(), pair.getRight()));
        }
        return compressedBins;
    }

    private static SimpleInterval makeCompressedBin(final List<SimpleInterval> bins, final int startIndex, final int endIndex) {
        if (startIndex < 0 || startIndex >= bins.size()) {
            throw new IllegalArgumentException("Invalid start bin index: " + startIndex);
        }
        if (endIndex < 0 || endIndex >= bins.size()) {
            throw new IllegalArgumentException("Invalid end bin index: " + endIndex);
        }
        final SimpleInterval startBin = bins.get(startIndex);
        final SimpleInterval endBin = bins.get(endIndex);
        return new SimpleInterval(startBin.getContig(), startBin.getStart(), endBin.getEnd());
    }

    /**
     * Compression function equivalent to the R function
     *
     * @param vals    array of values to compress
     * @param indices indices bounding compression intervals
     * @return compressed array as double[]
     */
    private static double[] compressCounts(final double[] vals, final List<Pair<Integer, Integer>> indices) {
        final double[] compressed = new double[indices.size()];
        int i = 0;
        for (final Pair<Integer, Integer> pair : indices) {
            final int startIndex = pair.getLeft();
            final int endIndex = pair.getRight();
            double sum = 0;
            for (int j = startIndex; j <= endIndex && j < vals.length; j++) {
                sum += vals[j];
            }
            compressed[i++] = sum;
        }
        return compressed;
    }

    protected static DepthMatrix normalizeMatrix(final DepthMatrix matrix, final Map<String, Double> sampleMedians) {
        if (matrix.isEmpty()) {
            return matrix;
        }
        final Map<String, double[]> normalizedMatrix = new HashMap<>();
        for (final String sample : matrix.getSampleSet()) {
            final double[] values = matrix.getSample(sample);
            if (!sampleMedians.containsKey(sample)) {
                throw new IllegalArgumentException("Sample " + sample + " has no median");
            }
            final double median = sampleMedians.get(sample);
            final double[] newValues = new double[values.length];
            for (int j = 0; j < matrix.getNumBins(); j++) {
                newValues[j] = values[j] / median;
            }
            normalizedMatrix.put(sample, newValues);
        }
        return new DepthMatrix(matrix.getBins(), normalizedMatrix);
    }
}
