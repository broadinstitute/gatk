package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

public final class DepthMatrix {

    private static final Median MEDIAN = new Median();

    private final List<SimpleInterval> bins;
    private final Map<String, double[]> matrix;

    public DepthMatrix(final List<SimpleInterval> bins, final Map<String, double[]> matrix) {
        Utils.nonNull(bins);
        Utils.nonNull(matrix);
        for (final double[] value : matrix.values()) {
            Utils.nonNull(value);
            Utils.validateArg(value.length == bins.size(), "Each vector must have same length as bins list");
        }
        this.bins = bins;
        this.matrix = matrix;
    }

    private static DepthMatrix fromListMap(final List<SimpleInterval> bins,
                                           final Map<String, List<Double>> coverageData) {
        final Map<String, double[]> arrayData = new HashMap<>();
        for (String sample : coverageData.keySet()) {
            final double[] values = coverageData.get(sample).stream().mapToDouble(Double::doubleValue).toArray();
            arrayData.put(sample, values);
        }
        return new DepthMatrix(bins, arrayData);
    }

    public static DepthMatrix fromDataSource(final SimpleInterval interval,
                                             final FeatureDataSource<DepthEvidence> dataSource) {
        final DepthMatrix depthMatrix = queryDataSource(dataSource, interval);

        // Remove bins that start or end exactly at query boundaries (mimics R behavior)
        final List<SimpleInterval> filteredBins = new ArrayList<>();
        final Map<String, List<Double>> filteredData = new HashMap<>();
        for (String sample : depthMatrix.getSampleSet()) {
            filteredData.put(sample, new ArrayList<>());
        }

        for (int i = 0; i < depthMatrix.getBins().size(); i++) {
            final SimpleInterval bin = depthMatrix.getBins().get(i);
            if (bin.getStart() <= interval.getEnd() && bin.getEnd() != interval.getStart()) {
                filteredBins.add(bin);
                for (String sample : depthMatrix.getSampleSet()) {
                    filteredData.get(sample)
                            .add(depthMatrix.getSample(sample)[i]);
                }
            }
        }

        // Fill gaps with zero-count bins
        return fillGapsInCoverageMatrix(filteredBins, filteredData);
    }

    private static DepthMatrix fillGapsInCoverageMatrix(final List<SimpleInterval> bins,
                                                        final Map<String, List<Double>> coverageData) {
        if (bins.size() <= 1) {
            return DepthMatrix.fromListMap(bins, coverageData);
        }

        final List<SimpleInterval> filledBins = new ArrayList<>(bins);
        final Map<String, List<Double>> filledData = new HashMap<>();

        // Initialize with existing data
        for (String sample : coverageData.keySet()) {
            filledData.put(sample, new ArrayList<>(coverageData.get(sample)));
        }

        final int binSize = bins.get(0).size(); // Assume uniform bin size

        for (int i = 0; i < bins.size() - 1; i++) {
            final SimpleInterval currentBin = bins.get(i);
            final SimpleInterval nextBin = bins.get(i + 1);

            final int gapLength = nextBin.getStart() - currentBin.getEnd() - 1;
            if (gapLength > 0) {
                // Fill gap with zero-count bins
                final String chr = currentBin.getContig();
                for (int pos = currentBin.getEnd(); pos < nextBin.getStart(); pos += binSize) {
                    final int gapEnd = Math.min(pos + binSize, nextBin.getStart());
                    final SimpleInterval gapBin = new SimpleInterval(chr, pos, gapEnd);

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
        return fromListMap(filledBins, filledData);
    }

    public static DepthMatrix fromDataSourceSubsampled(final SimpleInterval interval,
                                                       final FeatureDataSource<DepthEvidence> dataSource,
                                                       final int window,
                                                       final int numPoints,
                                                       final SAMSequenceDictionary dictionary) {
        final int padding = window / 2;
        final int pointSpacing = interval.getLengthOnReference() / numPoints;

        final List<SimpleInterval> pointIntervals = new ArrayList<>();
        for (int i = 0; i < numPoints; i++) {
            final int pointCenter = interval.getStart() + padding + (i * pointSpacing);
            final SimpleInterval pointInterval = new SimpleInterval(interval.getContig(), pointCenter, pointCenter);
            final SimpleInterval paddedInterval = pointInterval.expandWithinContig(padding, dictionary);
            pointIntervals.add(paddedInterval);
        }

        final List<SimpleInterval> bins = new ArrayList<>();
        final Map<String, List<Double>> listMatrix = new HashMap<>();
        for (int i = 0; i < numPoints; i++) {
            final SimpleInterval pointInterval = pointIntervals.get(i);
            final DepthMatrix pointMatrix = queryDataSource(dataSource, pointInterval);

            // Calculate median coverage for this window
            boolean hasData = false;
            for (String sample : pointMatrix.getSampleSet()) {
                final double[] values = pointMatrix.getSample(sample);
                if (!listMatrix.containsKey(sample)) {
                    listMatrix.put(sample, new ArrayList<>());
                }
                if (values.length > 0) {
                    listMatrix.get(sample).add(MEDIAN.evaluate(values));
                    hasData = true;
                }
            }
            if (hasData) {
                bins.add(pointInterval);
            }
        }
        return fromListMap(bins, listMatrix);
    }

    private static DepthMatrix queryDataSource(final FeatureDataSource<DepthEvidence> dataSource, final SimpleInterval interval) {
        final List<DepthEvidence> data = Lists.newArrayList(dataSource.query(interval));
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
            ++i;
        }
        return new DepthMatrix(bins, coverageData);
    }

    public double[] slice(final int binIndex, final Collection<String> samples) {
        Utils.validateArg(matrix.keySet().containsAll(samples), "Matrix does not contain all queried samples");
        if (matrix.isEmpty()) {
            return new double[0];
        }
        final double[] slice = new double[samples.size()];
        int j = 0;
        for (final String sample : samples) {
            slice[j++] = getSample(sample)[binIndex];
        }
        return slice;
    }

    public List<SimpleInterval> getBins() {
        return bins;
    }

    public double[] getSample(final String sample) {
        Utils.validateArg(matrix.containsKey(sample), "Matrix does not contain queried sample");
        return matrix.get(sample);
    }

    public boolean isEmpty() {
        return matrix.isEmpty();
    }

    public Set<String> getSampleSet() {
        return matrix.keySet();
    }

    public int getNumBins() {
        return bins.size();
    }
}
