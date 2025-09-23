package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.Utils;
import picard.sam.util.Pair;

import java.util.*;
import java.util.stream.Collectors;

public class DepthMatrixLoader {

    protected final FeatureDataSource<DepthEvidence> dataSource;
    protected final SAMSequenceDictionary dictionary;
    protected final int bins;
    protected final long largeSizeCutoff;
    protected final int largeVariantRegionPoints;
    protected final int largeVariantWindow;

    private record CompressionResult(double compression, SimpleInterval adjustedRegion) {
    }

    public DepthMatrixLoader(final FeatureDataSource<DepthEvidence> dataSource,
                             int bins, long largeSizeCutoff, int largeVariantRegionPoints, int largeVariantWindow,
                             final SAMSequenceDictionary dictionary) {
        Utils.nonNull(dataSource);
        Utils.validateArg(bins > 0, "bins must be greater than zero");
        Utils.validateArg(largeSizeCutoff >= 0, "large size cutoff must be non-negative");
        Utils.validateArg(largeVariantRegionPoints > 0, "large variant region points must be greater than zero");
        Utils.validateArg(largeVariantWindow >= 0, "large variant region window must be non-negative");
        this.dataSource = dataSource;
        this.dictionary = dictionary;
        this.bins = bins;
        this.largeSizeCutoff = largeSizeCutoff;
        this.largeVariantRegionPoints = largeVariantRegionPoints;
        this.largeVariantWindow = largeVariantWindow;
    }

    public DepthMatrix load(final SimpleInterval interval,
                            final Map<String, Double> sampleMedians) {
        Utils.nonNull(interval);
        if (!IntervalUtils.intervalIsOnDictionaryContig(interval, dictionary)) {
            throw new GATKException("Load interval " + interval + " failed to validate against the sequence dictionary");
        }
        Utils.nonNull(dataSource);
        Utils.nonNull(sampleMedians);

        final DepthMatrix rawData;
        if (interval.getLengthOnReference() > largeSizeCutoff) {
            rawData = loadLargeVariantMatrix(interval, dataSource);
        } else {
            rawData = loadStandardVariantMatrix(interval, dataSource);
        }

        // Calculate compression and region bounds
        final CompressionResult compressionResult = calculateCompression(rawData, interval, bins);

        // Trim bins from compression
        final DepthMatrix trimmedData = trimBins(rawData, compressionResult);

        // Set 0 entries to 1
        setZeroesToOnes(trimmedData);

        final DepthMatrix compressedData;
        final Map<String, Double> compressedMedians;
        if (compressionResult == null || compressionResult.compression <= 1) {
            // No compression needed
            compressedData = trimmedData;
            compressedMedians = sampleMedians;
        } else {
            // Apply compression and normalization
            compressedData = compressMatrix(trimmedData, compressionResult);
            compressedMedians = compressMedians(sampleMedians, compressionResult);
        }
        final DepthMatrix normalizedData = normalizeData(compressedData, compressedMedians);
        return trimMatrix(normalizedData);
    }

    private static void setZeroesToOnes(final DepthMatrix matrix) {
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

    private DepthMatrix loadLargeVariantMatrix(final SimpleInterval interval,
                                               final FeatureDataSource<DepthEvidence> dataSource) {
        int padding = largeVariantWindow / 2;
        int pointSpacing = interval.getLengthOnReference() / largeVariantRegionPoints;

        List<SimpleInterval> pointIntervals = new ArrayList<>();
        for (int i = 0; i < largeVariantRegionPoints; i++) {
            int pointCenter = interval.getStart() + padding + (i * pointSpacing);
            int pointStart = Math.max(1, pointCenter - padding);
            int pointEnd = Math.min(pointCenter + padding + 1, dictionary.getSequence(interval.getContig()).getSequenceLength());
            final SimpleInterval newInterval = new SimpleInterval(interval.getContig(), pointStart, pointEnd);
            if (!IntervalUtils.intervalIsOnDictionaryContig(newInterval, dictionary)) {
                throw new GATKException("Padded interval " + newInterval + " failed to validate against the sequence dictionary");
            }
            pointIntervals.add(newInterval);
        }

        final List<SimpleInterval> bins = new ArrayList<>();
        final Map<String, List<Double>> listMatrix = new HashMap<>();

        final Median median = new Median();
        for (int i = 0; i < largeVariantRegionPoints; i++) {
            final SimpleInterval pointInterval = pointIntervals.get(i);
            final DepthMatrix windowCov = queryDataSource(dataSource, pointInterval);

            // Calculate median coverage for this window
            boolean hasData = false;
            for (String sample : windowCov.getSampleSet()) {
                double[] values = windowCov.getSample(sample);
                if (!listMatrix.containsKey(sample)) {
                    listMatrix.put(sample, new ArrayList<>());
                }
                if (values.length > 0) {
                    listMatrix.get(sample).add(median.evaluate(values));
                    hasData = true;
                }
            }
            if (hasData) {
                bins.add(pointInterval);
            }
        }
        final Map<String, double[]> arrayMatrix = new HashMap<>();
        for (String sample : listMatrix.keySet()) {
            arrayMatrix.put(sample, listMatrix.get(sample).stream().mapToDouble(Double::doubleValue).toArray());
        }
        return new DepthMatrix(bins, arrayMatrix);
    }

    private DepthMatrix loadStandardVariantMatrix(final SimpleInterval interval,
                                                  final FeatureDataSource<DepthEvidence> dataSource) {
        DepthMatrix depthMatrix = queryDataSource(dataSource, interval);

        // Remove bins that start or end exactly at query boundaries (mimics R behavior)
        List<SimpleInterval> filteredBins = new ArrayList<>();
        Map<String, List<Double>> filteredData = new HashMap<>();
        for (String sample : depthMatrix.getSampleSet()) {
            filteredData.put(sample, new ArrayList<>());
        }

        for (int i = 0; i < depthMatrix.getBins().size(); i++) {
            SimpleInterval bin = depthMatrix.getBins().get(i);
            if (bin.getStart() <= interval.getEnd() && bin.getEnd() != interval.getStart()) {
                filteredBins.add(bin);
                for (String sample : depthMatrix.getSampleSet()) {
                    filteredData.get(sample)
                            .add(depthMatrix.getSample(sample)[i]);
                }
            }
        }

        // Fill gaps with zero-count bins
        return fillGapsInCoverageMatrix(filteredBins, filteredData, interval.getContig());
    }

    private DepthMatrix fillGapsInCoverageMatrix(final List<SimpleInterval> bins,
                                                 final Map<String, List<Double>> coverageData,
                                                 final String chr) {
        if (bins.size() <= 1) {
            return createCoverageMatrixFromLists(bins, coverageData);
        }

        List<SimpleInterval> filledBins = new ArrayList<>(bins);
        Map<String, List<Double>> filledData = new HashMap<>();

        // Initialize with existing data
        for (String sample : coverageData.keySet()) {
            filledData.put(sample, new ArrayList<>(coverageData.get(sample)));
        }

        int binSize = bins.get(0).size(); // Assume uniform bin size

        for (int i = 0; i < bins.size() - 1; i++) {
            SimpleInterval currentBin = bins.get(i);
            SimpleInterval nextBin = bins.get(i + 1);

            int gapLength = nextBin.getStart() - currentBin.getEnd() - 1;
            if (gapLength > 0) {
                // Fill gap with zero-count bins
                for (int pos = currentBin.getEnd(); pos < nextBin.getStart(); pos += binSize) {
                    int gapEnd = Math.min(pos + binSize, nextBin.getStart());
                    SimpleInterval gapBin = new SimpleInterval(chr, pos, gapEnd);

                    int insertIndex = filledBins.size();
                    for (int j = 0; j < filledBins.size(); j++) {
                        if (filledBins.get(j).getStart() > pos) {
                            insertIndex = j;
                            break;
                        }
                    }

                    filledBins.add(insertIndex, gapBin);

                    // Add zero values for all samples
                    for (String sample : filledData.keySet()) {
                        filledData.get(sample).add(insertIndex, 0.0);
                    }
                }
            }
        }

        return createCoverageMatrixFromLists(filledBins, filledData);
    }

    private DepthMatrix trimBins(DepthMatrix cov1, final CompressionResult compressionResult) {
        if (compressionResult == null) {
            return cov1;
        }
        final int trimStart = compressionResult.adjustedRegion.getStart();
        final int trimEnd = compressionResult.adjustedRegion.getEnd();
        Set<Integer> indicesToRemove = new HashSet<>();
        for (int i = 0; i < cov1.getBins().size(); i++) {
            final SimpleInterval bin = cov1.getBins().get(i);
            if (bin.getEnd() <= trimStart || bin.getStart() >= trimEnd) {
                indicesToRemove.add(i);
            }
        }
        if (indicesToRemove.isEmpty()) {
            return cov1;
        }

        final int newNumBins = cov1.getBins().size() - indicesToRemove.size();
        final List<SimpleInterval> newBins = new ArrayList<>(newNumBins);
        for (int i = 0; i < cov1.getBins().size(); i++) {
            if (!indicesToRemove.contains(i)) {
                newBins.add(cov1.getBins().get(i));
            }
        }
        final Map<String, double[]> newCoverageData = new HashMap<>();
        for (final String sample : cov1.getSampleSet()) {
            final double[] values = cov1.getSample(sample);
            final double[] newValues = new double[newNumBins];
            int j = 0;
            for (int i = 0; i < values.length; i++) {
                if (!indicesToRemove.contains(i)) {
                    newValues[j++] = values[i];
                }
            }
            newCoverageData.put(sample, newValues);
        }
        return new DepthMatrix(newBins, newCoverageData);
    }

    private DepthMatrix trimMatrix(final DepthMatrix matrix) {
        final int numBins = matrix.getNumBins();
        final boolean subset = numBins >= 4;
        final List<SimpleInterval> newBins = subset ? matrix.getBins().subList(1, numBins - 1) : matrix.getBins();
        final Map<String, double[]> newMatrix = new HashMap<>();
        for (final String sample : matrix.getSampleSet()) {
            final double[] counts = matrix.getSample(sample);
            final double[] trimmedCounts;
            if (subset) {
                trimmedCounts = Arrays.copyOfRange(counts, 1, numBins - 1);
            } else {
                trimmedCounts = counts;
            }
            newMatrix.put(sample, trimmedCounts);
        }
        return new DepthMatrix(newBins, newMatrix);
    }

    private static CompressionResult calculateCompression(DepthMatrix cov1, SimpleInterval queryInterval, int bins) {
        int start = queryInterval.getStart() - 1; // covert to 0-based to match R code
        int end = queryInterval.getEnd();

        // Find bins that overlap with the query interval
        int startBinIdx = -1, endBinIdx = -1;

        for (int i = 0; i < cov1.getBins().size(); i++) {
            SimpleInterval bin = cov1.getBins().get(i);
            if (startBinIdx == -1 && bin.getStart() - 1 >= start) {
                startBinIdx = i;
            }
            if (bin.getEnd() <= end) {
                endBinIdx = i;
            }
        }

        if (startBinIdx == -1) startBinIdx = 0;
        if (endBinIdx == -1) endBinIdx = cov1.getBins().size() - 1;

        int numInternalBins = endBinIdx - startBinIdx + 1;

        if (numInternalBins < bins) {
            bins = numInternalBins;
            if (bins == 0) {
                return null;
            } else if (bins <= 1) {
                SimpleInterval adjustedRegion = new SimpleInterval(
                        queryInterval.getContig(),
                        cov1.getBins().get(0).getStart(),
                        cov1.getBins().get(cov1.getBins().size() - 1).getEnd()
                );
                return new CompressionResult(1, adjustedRegion);
            }
        }

        // Adjust for even compression
        int remainderForRemoval = (int) Math.floor(((double) (numInternalBins % bins)) / 2.0);
        int newStartBinIdx = startBinIdx + remainderForRemoval;
        int newEndBinIdx = endBinIdx - remainderForRemoval;

        // Calculate compression factor
        double compression = (newEndBinIdx - newStartBinIdx + 1) / (double) bins;

        SimpleInterval adjustedRegion = new SimpleInterval(
                queryInterval.getContig(),
                cov1.getBins().get(newStartBinIdx).getStart(),
                cov1.getBins().get(newEndBinIdx).getEnd()
        );

        return new CompressionResult(compression, adjustedRegion);
    }

    private DepthMatrix compressMatrix(DepthMatrix cov1,
                                       CompressionResult compressionResult) {
        if (cov1.isEmpty()) {
            return cov1;
        }
        final String firstSample = cov1.getSampleSet().iterator().next();
        final List<Pair<Integer, Integer>> compressionIndices = computeCompressionIndices(cov1.getSample(firstSample), compressionResult.compression);
        final List<SimpleInterval> compressedBins = compressBins(cov1.getBins(), compressionIndices);
        Map<String, double[]> compressedMatrix = new LinkedHashMap<>();
        List<String> sampleNames = new ArrayList<>(cov1.getSampleSet());
        for (final String sample : sampleNames) {
            double[] columnValues = cov1.getSample(sample);

            // Apply compression function to this column
            double[] compressedColumn = compressColumn(columnValues, compressionIndices);
            compressedMatrix.put(sample, compressedColumn);
        }
        return new DepthMatrix(compressedBins, compressedMatrix);
    }

    private Map<String, Double> compressMedians(Map<String, Double> allnorm,
                                                CompressionResult compressionResult) {
        // Approximate rebinned per-sample medians
        Map<String, Double> adjustedMedians = new HashMap<>();
        for (String sample : allnorm.keySet()) {
            double median = allnorm.get(sample);
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
    private static List<Pair<Integer, Integer>> computeCompressionIndices(double[] vals, double compression) {
        Utils.validateArg(compression > 1., "Compression must be greater than 1");
        if (vals.length == 0) return Collections.emptyList();
        List<Pair<Integer, Integer>> result = new ArrayList<>();

        // Convert compression to int for array indexing
        final int firstEnd = Math.min(((int) compression) - 1, vals.length - 1);
        result.add(new Pair<>(0, firstEnd));

        // Process remaining chunks
        int numChunks = (int) (vals.length / compression);
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

    private List<SimpleInterval> compressBins(List<SimpleInterval> bins, List<Pair<Integer, Integer>> indices) {
        final List<SimpleInterval> compressedBins = new ArrayList<>();
        for (final Pair<Integer, Integer> pair : indices) {
            final int startIndex = pair.getLeft();
            final int endIndex = pair.getRight();
            final SimpleInterval startBin = bins.get(startIndex);
            final SimpleInterval endBin = bins.get(endIndex);
            final SimpleInterval compressedBin = new SimpleInterval(startBin.getContig(), startBin.getStart(), endBin.getEnd());
            if (!IntervalUtils.intervalIsOnDictionaryContig(compressedBin, dictionary)) {
                throw new GATKException("Compressed interval " + compressedBin + " failed to validate against the sequence dictionary");
            }
            compressedBins.add(compressedBin);
        }
        return compressedBins;
    }

    /**
     * Compression function equivalent to the R function
     *
     * @param vals    array of values to compress
     * @param indices indices bounding compression intervals
     * @return compressed array as double[]
     */
    private static double[] compressColumn(double[] vals, List<Pair<Integer, Integer>> indices) {
        double[] newcol = new double[indices.size()];
        int i = 0;
        for (final Pair<Integer, Integer> pair : indices) {
            final int startIndex = pair.getLeft();
            final int endIndex = pair.getRight();
            double sum = 0;
            for (int j = startIndex; j <= endIndex && j < vals.length; j++) {
                sum += vals[j];
            }
            newcol[i++] = sum;
        }
        return newcol;
    }


    private static DepthMatrix normalizeData(DepthMatrix matrix, Map<String, Double> medians) {
        List<String> samples = new ArrayList<>(matrix.getSampleSet());
        if (samples.isEmpty()) {
            return matrix;
        }
        int numBins = matrix.getSample(samples.get(0)).length;
        Map<String, double[]> normalizedMatrix = new HashMap<>(SVUtils.hashMapCapacity(matrix.getSampleSet().size()));
        for (final String sample : samples) {
            double[] values = matrix.getSample(sample);
            double median = medians.getOrDefault(sample, 1.0);
            final double[] newValues = new double[values.length];
            for (int j = 0; j < numBins; j++) {
                newValues[j] = values[j] / median;
            }
            normalizedMatrix.put(sample, newValues);
        }
        return new DepthMatrix(matrix.getBins(), normalizedMatrix);
    }

    // Helper methods
    private static DepthMatrix queryDataSource(FeatureDataSource<DepthEvidence> dataSource, SimpleInterval interval) {
        final Iterator<DepthEvidence> iter = dataSource.query(interval);
        final List<DepthEvidence> data = Lists.newArrayList(iter);
        final List<SimpleInterval> bins = data.stream().map(d -> new SimpleInterval(d.getContig(), d.getStart(), d.getEnd())).collect(Collectors.toList());
        final SVFeaturesHeader header = (SVFeaturesHeader) dataSource.getHeader();
        final List<String> sampleNames = header.getSampleNames();
        final Map<String, double[]> coverageData = new HashMap<>();
        for (final String sampleName : sampleNames) {
            coverageData.put(sampleName, new double[bins.size()]);
        }
        int i = 0;
        for (final DepthEvidence depthEvidence : data) {
            int j = 0;
            for (final String sampleName : sampleNames) {
                coverageData.get(sampleName)[i] = depthEvidence.getCounts()[j++];
            }
            i += 1;
        }

        // This would be replaced with actual GATK FeatureDataSource query logic
        return new DepthMatrix(bins, coverageData);
    }

    private static DepthMatrix createCoverageMatrixFromLists(List<SimpleInterval> bins,
                                                             Map<String, List<Double>> coverageData) {
        Map<String, double[]> arrayData = new HashMap<>();
        for (String sample : coverageData.keySet()) {
            List<Double> values = coverageData.get(sample);
            arrayData.put(sample, values.stream().mapToDouble(Double::doubleValue).toArray());
        }
        return new DepthMatrix(bins, arrayData);
    }
}
