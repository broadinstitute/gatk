package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.solver.SimulatedAnnealingSolver;
import scala.Tuple2;

import java.util.*;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;

public final class ReadDepthModel {

    private final List<ReadDepthCluster> clusteredEvents;
    private final NormalDistribution standardNormal;
    private ReadDepthModelParameters parameters;
    public ReadDepthModel(final SVIntervalTree<LargeSimpleSV> eventsTree, final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        this.parameters = new ReadDepthModelParameters();
        this.standardNormal = new NormalDistribution(0, 1);
        setSamplerSeed(0);
        this.clusteredEvents = clusterEvents(eventsTree, copyRatioSegmentOverlapDetector, dictionary);
    }

    private static double unscaledNormal(final double x, final double mu, final double sigma) {
        final double diff = (x - mu) / sigma;
        return Math.exp(-0.5 * diff * diff) / sigma;
    }

    private static double unscaledLogNormal(final double x, final double mu, final double sigma) {
        final double diff = (x - mu) / sigma;
        return -0.5 * diff * diff - Math.log(sigma);
    }

    public void setParameters(final ReadDepthModelParameters parameters) {
        this.parameters = parameters;
    }

    public void setSamplerSeed(final long seed) {
        standardNormal.reseedRandomGenerator(seed);
    }

    private List<ReadDepthCluster> clusterEvents(final SVIntervalTree<LargeSimpleSV> eventsTree, final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        final Set<LargeSimpleSV> visited = new HashSet<>(SVUtils.hashMapCapacity(eventsTree.size()));
        final List<List<ReadDepthEvent>> clusteredEvents = new ArrayList<>(eventsTree.size());
        final Iterator<SVIntervalTree.Entry<LargeSimpleSV>> iter = eventsTree.iterator();
        final Set<LargeSimpleSV> cluster = new HashSet<>(SVUtils.hashMapCapacity(eventsTree.size()));
        while (iter.hasNext()) {
            final SVIntervalTree.Entry<LargeSimpleSV> entry = iter.next();
            final LargeSimpleSV event = entry.getValue();
            if (visited.contains(event)) continue;
            cluster.clear();
            cluster.add(event);
            int oldClusterSize;
            do {
                final List<LargeSimpleSV> oldCluster = new ArrayList<>(cluster); //Have to copy the set to avoid ConcurrentModificationException
                oldClusterSize = cluster.size();
                for (final LargeSimpleSV clusterEvent : oldCluster) {
                    cluster.addAll(Utils.stream(eventsTree.overlappers(clusterEvent.getInterval())).map(SVIntervalTree.Entry::getValue).collect(Collectors.toList()));
                }
            } while (cluster.size() > oldClusterSize);
            visited.addAll(cluster);

            final List<ReadDepthEvent> newCluster = new ArrayList<>(cluster.size());
            int eventId = 0;
            final Iterator<LargeSimpleSV> newClusterIter = cluster.iterator();
            while (newClusterIter.hasNext()) {
                newCluster.add(new ReadDepthEvent(eventId, newClusterIter.next()));
                eventId++;
            }
            clusteredEvents.add(newCluster);
        }
        return clusteredEvents.stream().map(events -> new ReadDepthCluster(events, copyRatioSegmentOverlapDetector, dictionary))
                .collect(Collectors.toList());
    }

    public Collection<ReadDepthEvent> solve() {
        return clusteredEvents.stream().flatMap(cluster -> solve(cluster).stream()).collect(Collectors.toList());
    }

    private Collection<ReadDepthEvent> solve(final ReadDepthCluster cluster) {
        final List<ReadDepthEvent> events = cluster.getEventsList();
        final Function<double[], Double> energyFunction = x -> -computeLogPosterior(x, cluster);
        final int size = events.size();
        final Supplier<double[]> sampler = () -> standardNormal.sample(size);
        final double[] lowerBound = new double[size];
        final double[] upperBound = new double[size];
        Arrays.fill(lowerBound, 0);
        Arrays.fill(upperBound, parameters.maxPloidy);
        final SimulatedAnnealingSolver sa = new SimulatedAnnealingSolver(size, energyFunction, sampler, lowerBound, upperBound);

        final double[] x0 = new double[size];
        final double T0 = parameters.saT0;
        final int numSteps = parameters.saSteps;
        final double finalEnergy = sa.solve(x0, T0, numSteps);

        final double[] minEnergyState = sa.getSolution();
        for (int i = 0; i < events.size(); i++) {
            events.get(i).setState(minEnergyState[i]);
        }
        return events;

        /*final GradientDescentSolver g2 = new GradientDescentSolver(parameters, eventsList.size());
        g2.initialize(sa.r, sa.q);

        int numIter = 0;
        int numDeltaZero = 0;
        double deltaRiTotal = 0;
        boolean lastDeltaZero = false;
        while (numDeltaZero < 100 && numIter < parameters.maxIter) {
            if (numIter % 10000 == 0) {
                logger.info("\titer = " + numIter);
                final StringBuilder stringBuilder = new StringBuilder();
                for (int i = 0; i < g2.r.length; i++) {
                    stringBuilder.append("(");
                    stringBuilder.append(g2.r[i]);
                    stringBuilder.append(", ");
                    stringBuilder.append(g2.q[i]);
                    stringBuilder.append(") ");
                }
                logger.info("(ri, qi) = {" + stringBuilder.toString() + "}");
            }
            deltaRiTotal = g2.computeGradient(eventsList, nearestCallDistances, copyNumberInfo);
            g2.computeDelta(numIter);
            g2.step();

            if (deltaRiTotal < parameters.absoluteTolerance) {
                if (lastDeltaZero) {
                    numDeltaZero++;
                }
                lastDeltaZero = true;
            } else {
                numDeltaZero = 0;
                lastDeltaZero = false;
            }
            numIter++;
        }
        */
        /*if (deltaRiTotal > parameters.absoluteTolerance) {
            final double qtemp = g2.q[0];
            for (int i = -100; i <= 2100; i++) {
                g2.q[0] = i / 1000.0;
                final double energy = computeLogLikelihood(parameters, g2.r, g2.q, copyNumberInfo, eventsList, nearestCallDistances);
                System.out.println(g2.q[0] + "\t" + energy);
            }
            g2.q[0] = qtemp;
        }*/
        /*
        final double logPosterior = g2.computeLogPosterior(eventsList, nearestCallDistances, copyNumberInfo);
        logger.info("Finished after " + numIter + " iterations with log posterior " + logPosterior + " and delta " + deltaRiTotal);
        final StringBuilder stringBuilder = new StringBuilder();
        for (int i = 0; i < g2.r.length; i++) {
            stringBuilder.append("(");
            stringBuilder.append(g2.r[i]);
            stringBuilder.append(", ");
            stringBuilder.append(g2.q[i]);
            stringBuilder.append(") ");
        }
        logger.info("(ri, qi) = {" + stringBuilder.toString() + "}");
        final List<ReadDepthEvent> finalCalls = new ArrayList<>(eventsList);
        for (final ReadDepthEvent call : finalCalls) {
            call.setResults(g2.r[call.id], g2.q[call.id]);
        }
        return finalCalls;
        */
    }

    private double computeLogPosterior(final double[] x, final ReadDepthCluster cluster) {
        return computeLogLikelihood(x, cluster) + computeLogPrior(x, cluster);
    }

    private double computeLogLikelihood(final double[] x, final ReadDepthCluster cluster) {
        return cluster.getCopyNumberInfo().stream().mapToDouble(entry -> computeCopyNumberLikelihood(x, entry)).sum()
                + cluster.getEventsList().stream().mapToDouble(event -> computeReadEvidenceLikelihood(x, event)).sum()
                + computeCallDistanceLikelihood(x, cluster);
    }

    private double computeLogPrior(final double[] x, final ReadDepthCluster cluster) {
        return cluster.getCopyNumberInfo().stream().mapToDouble(entry -> computePloidyPrior(x, entry)).sum();
    }

    private double computeCallDistanceLikelihood(final double[] x, final ReadDepthCluster cluster) {
        double total = 0;
        final List<ReadDepthEvent> events = cluster.getEventsList();
        final List<Tuple2<Integer, Integer>> nearestCallDistances = cluster.getNearestCallDistances();
        for (int i = 0; i < events.size(); i++) {
            final double std = parameters.meanInsertSize / Math.max(parameters.callDistancePseudocount, parameters.callDistancePseudocount + stateToProbability(x[i]));
            total += events.get(i).getEvent().getSize() * (unscaledLogNormal(nearestCallDistances.get(i)._1, 0, std) + unscaledLogNormal(nearestCallDistances.get(i)._2, 0, std));
        }
        return total;
    }

    private double computeReadEvidenceLikelihood(final double[] x, final ReadDepthEvent event) {
        final double expectedEvidence = x[event.getId()] * parameters.copyNeutralDepth * parameters.expectedReadEvidenceFraction;
        final double sigma = parameters.expectedReadEvidenceStd * parameters.copyNeutralDepth;
        return event.getEvent().getSize() * (unscaledLogNormal(event.getEvent().getReadPairEvidence(), expectedEvidence, sigma)
                + unscaledLogNormal(event.getEvent().getSplitReadEvidence(), expectedEvidence, sigma));
    }

    private double computeCopyNumberLikelihood(final double[] x, final Tuple2<List<OverlapInfo>, Double> entry) {
        final List<OverlapInfo> overlapInfoList = entry._1;
        final double calledCopyNumber = entry._2;
        return overlapInfoList.stream()
                .mapToDouble(info -> {
                    final double estimatedCopyNumber = 2 + info.idsAndCoefficients.stream().mapToDouble(tuple -> tuple._2 * x[tuple._1]).sum();
                    return unscaledLogNormal(calledCopyNumber, estimatedCopyNumber, parameters.copyNumberStd) * info.size;
                }).sum();
    }

    private double stateToProbability(final double q) {
        return 1.0 - Math.exp(-q / 0.5);
    }

    private double computePloidyPrior(final double[] x, final Tuple2<List<OverlapInfo>, Double> entry) { //, final double[] deltaR) {
        final List<OverlapInfo> overlapInfoList = entry._1;
        double total = 0;
        for (final OverlapInfo info : overlapInfoList) {
            double totalPloidy = info.idsAndCoefficients.stream().mapToDouble(pair -> x[pair._1]).sum();
            final double estimatedPloidy;
            if (totalPloidy > parameters.maxPloidy) {
                estimatedPloidy = parameters.maxPloidy - totalPloidy;
            } else if (totalPloidy < 0) {
                estimatedPloidy = totalPloidy;
            } else {
                estimatedPloidy = 0;
            }
            total += unscaledLogNormal(estimatedPloidy, 0, parameters.parameterConstraintStd) * info.size;
        }
        return total;
    }

    private final static class ReadDepthModelParameters {
        /*
        public final double gradientDelta = 1e-4;
        public final double learningRate = 1e-3; //1e-9;
        public final double absoluteTolerance = 1e-4;
        public final int maxIter = 200000;
        public final double maxStepSize = 0.001;

        public final double maxR = 1;
        */

        public final int saSteps = 10000;
        public final double saT0 = 10000;

        public final double maxPloidy = 2;
        public final double copyNeutralDepth = 30;
        public final double meanInsertSize = 100;
        public final double parameterConstraintStd = 0.01;
        public final double pqDifferenceStd = 1;
        public final double expectedReadEvidenceFraction = 0.5;
        public final double expectedReadEvidenceStd = 0.1;
        public final double copyNumberStd = 0.1;
        public final double callDistancePseudocount = 1e-3;
    }

    static final class OverlapInfo {
        public List<Tuple2<Integer, Integer>> idsAndCoefficients;
        public double size;

        public OverlapInfo(List<Tuple2<Integer, Integer>> idsAndCoefficients, double size) {
            this.idsAndCoefficients = idsAndCoefficients;
            this.size = size;
        }
    }

    /*
    private double probabilityGivenZygosityPrior(final double[] x, final ReadDepthCluster cluster) {
        double total = 0;
        for (int i = 0; i < r.length; i++) {
            total += unscaledLogNormal(Math.min(1, q[i]) - r[i], 0, parameters.pqDifferenceStd) * events.get(i).getSize();
        }
        return total;
    }
    */

    /*
    private double quantizationPrior(final double[] x, final int numStates, final ReadDepthCluster cluster) {
        double total = 0;
        for (int i = 0; i < x.length; i++) {
            final double size = events.get(i).getSize();
            for (int j = 0; j < numStates; j++) {
                total += unscaledLogNormal(x[i], j, 1.0) * size;
            }
        }
        return total;
    }
    */

    /*
    private double parameterBoundsPrior(final double[] x, final double xMax, final ReadDepthCluster cluster) {
        double total = 0;
        for (int i = 0; i < x.length; i++) {
            final double y;
            if (x[i] > xMax) {
                y = xMax - x[i];
            } else if (x[i] < 0) {
                y = x[i];
            } else {
                y = 0;
            }
            total += unscaledLogNormal(y, 0, parameterConstraintStd) * events.get(i).getSize();
        }
        return total;
    }
    */
}
