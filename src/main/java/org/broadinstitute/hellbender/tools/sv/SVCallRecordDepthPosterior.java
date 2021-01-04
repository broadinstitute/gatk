package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.util.IntervalTree;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

public class SVCallRecordDepthPosterior {
    private final SVCallRecord record;
    private final int posteriorIntervalsOverlap;
    private final Map<String, CopyNumberPosteriorDistribution> samplePosteriors;

    public SVCallRecordDepthPosterior(final SVCallRecord record,
                                      final int posteriorIntervalsOverlap,
                                      final Map<String, CopyNumberPosteriorDistribution> samplePosteriors) {
        this.record = record;
        this.posteriorIntervalsOverlap = posteriorIntervalsOverlap;
        this.samplePosteriors = samplePosteriors;
    }

    public SVCallRecord getRecord() {
        return record;
    }

    public int getPosteriorIntervalsOverlap() {
        return posteriorIntervalsOverlap;
    }

    public Map<String, CopyNumberPosteriorDistribution> getSamplePosteriors() {
        return samplePosteriors;
    }

    public static SVCallRecordDepthPosterior create(final List<IntegerCopyNumberState> states,
                                                    final Collection<String> samples,
                                                    final SVCallRecord record,
                                                    final List<IntervalTree<Map<String,double[]>>> currentPosteriorsTreeList) {
        final SimpleInterval interval = new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB());
        final List<IntervalTree.Node<Map<String,double[]>>> posteriorsList = getBestPosteriors(interval, currentPosteriorsTreeList);
        final int overlap = getOverlap(interval, posteriorsList);
        final Map<String, CopyNumberPosteriorDistribution> samplePosteriors = samples.stream()
                .map(sample -> new HashMap.SimpleImmutableEntry<>(sample, getIntervalPosterior(states, record.getContigA(), interval, posteriorsList, sample)))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
        return new SVCallRecordDepthPosterior(record, overlap, samplePosteriors);
    }

    private static int getOverlap(final SimpleInterval variantInterval,
                                  final List<IntervalTree.Node<Map<String,double[]>>> posteriorsList) {
        final String contig = variantInterval.getContig();
        return posteriorsList.stream().map(node -> new SimpleInterval(contig, node.getStart(), node.getEnd()))
                .mapToInt(interval -> interval.intersect(variantInterval).size())
                .sum();
    }


    private static CopyNumberPosteriorDistribution getIntervalPosterior(final List<IntegerCopyNumberState> states,
                                                                        final String contig,
                                                                        final SimpleInterval variantInterval,
                                                                        final List<IntervalTree.Node<Map<String,double[]>>> posteriorsList,
                                                                        final String sample) {
        final int numCopyStates = states.size();
        final double[] copyStateSums = new double[numCopyStates];
        Arrays.fill(copyStateSums, Double.MIN_VALUE);
        int overlapSize = 0;
        for (final IntervalTree.Node<Map<String,double[]>> node : posteriorsList) {
            final SimpleInterval posteriorInterval = new SimpleInterval(contig, node.getStart(), node.getEnd());
            final int overlap = posteriorInterval.intersect(variantInterval).size();
            final double overlapFraction = overlap / (double) posteriorInterval.getLengthOnReference();
            overlapSize += overlap;
            final double[] dist = node.getValue().get(sample);
            for (int j = 0; j < numCopyStates; j++) {
                copyStateSums[j] += dist[j] * overlapFraction;
            }
        }
        if (overlapSize == 0) {
            return defaultDistribution(states);
        }

        double denom = 0;
        final double maxStateSum = DoubleStream.of(copyStateSums).max().getAsDouble();
        for (int i = 0; i < copyStateSums.length; i++) {
            // Normalize to avoid underflow error
            copyStateSums[i] -= maxStateSum;
            denom += FastMath.exp(copyStateSums[i]);
        }
        final double logDenom = Math.log(denom);
        final Map<IntegerCopyNumberState,Double> eventPosterior = new HashMap<>(SVUtils.hashMapCapacity(numCopyStates));
        for (int i = 0; i < copyStateSums.length; i++) {
            final Double p = copyStateSums[i] - logDenom;
            eventPosterior.put(states.get(i), p);
        }
        return new CopyNumberPosteriorDistribution(eventPosterior);
    }

    private static CopyNumberPosteriorDistribution defaultDistribution(final List<IntegerCopyNumberState> states) {
        final double p = FastMath.log(states.size() == 0 ? 1.0 : 1.0 / states.size());
        final Map<IntegerCopyNumberState, Double> logDist = states.stream()
                .map(state -> new AbstractMap.SimpleImmutableEntry<>(state, p))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
        return new CopyNumberPosteriorDistribution(logDist);
    }

    private static List<IntervalTree.Node<Map<String,double[]>>> getBestPosteriors(final SimpleInterval interval,
                                                                                   final List<IntervalTree<Map<String,double[]>>> currentPosteriorsTreeList) {
        int maxOverlap = -1;
        int maxIntervalCount = 0;
        List<IntervalTree.Node<Map<String,double[]>>> maxNodeList = null;
        for (int i = 0; i < currentPosteriorsTreeList.size(); i++) {
            final List<IntervalTree.Node<Map<String,double[]>>> nodeList = Lists.newArrayList(currentPosteriorsTreeList.get(i)
                    .overlappers(interval.getStart(), interval.getEnd()));
            final int overlap = nodeList.stream()
                    .map(node -> new SimpleInterval(interval.getContig(), node.getStart(), node.getEnd()))
                    .mapToInt(nodeInterval -> nodeInterval.intersect(interval).getLengthOnReference())
                    .sum();
            if (overlap > maxOverlap || (overlap == maxOverlap && nodeList.size() > maxIntervalCount)) {
                maxOverlap = overlap;
                maxIntervalCount = nodeList.size();
                maxNodeList = nodeList;
            }
        }
        return maxNodeList;
    }

}
