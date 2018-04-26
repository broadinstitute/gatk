package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.solver.GradientDescentSolver;
import org.broadinstitute.hellbender.utils.solver.SimulatedAnnealingSolver;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.TYPES.DEL;

public final class ReadDepthModel implements Serializable {

    public static final long serialVersionUID = 1L;

    public static final double MAX_PLOIDY = 2;
    public static final double COPY_NEUTRAL_DEPTH = 30;
    public static final int LINK_DENSITY_WIDTH = 50000;
    public static final int MAX_DUPLICATION_NUMBER = 5;

    private final Map<SimpleSVType.TYPES,List<ReadDepthCluster>> clusteredEvents;
    private ReadDepthModelParameters parameters;
    private final Logger logger = LogManager.getLogger(this.getClass());
    private long seed;

    public ReadDepthModel(final SVIntervalTree<LargeSimpleSV> callableEventsTree, final SVIntervalTree<EvidenceTargetLink> links, final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SVIntervalTree<Object> mappableIntervalTree, final SAMSequenceDictionary dictionary) {
        this.parameters = new ReadDepthModelParameters();
        setSamplerSeed(0);
        this.clusteredEvents = clusterEvents(callableEventsTree, copyRatioSegmentOverlapDetector, dictionary);
        setLinkDensity(links, dictionary);
        setMappabilityIndex(mappableIntervalTree);
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
        this.seed = seed;
    }

    public void setMappabilityIndex(final SVIntervalTree<Object> mappableTree) {
        final Collection<ReadDepthEvent> events = getEvents();
        for (final ReadDepthEvent event : events) {
            final SVInterval eventInterval = event.getEvent().getInterval();
            final int mappableOverlap = Utils.stream(mappableTree.overlappers(eventInterval)).map(SVIntervalTree.Entry::getInterval)
                    .mapToInt(interval -> interval.overlapLen(eventInterval)).sum();
            event.mappabilityIndex = mappableOverlap / (double) eventInterval.getLength();
        }
    }

    public void setLinkDensity(final SVIntervalTree<EvidenceTargetLink> links, final SAMSequenceDictionary dictionary) {
        final Collection<ReadDepthEvent> events = getEvents();
        for (final ReadDepthEvent event : events) {
            final SVInterval eventInterval = event.getEvent().getInterval();
            final Set<EvidenceTargetLink> eventLinks = new HashSet<>(event.getEvent().getSupportingEvidence());
            final SVInterval startPoint = new SVInterval(eventInterval.getContig(), eventInterval.getStart(), eventInterval.getStart());
            final SVInterval endPoint = new SVInterval(eventInterval.getContig(), eventInterval.getEnd(), eventInterval.getEnd());
            final SVInterval startInterval = SVIntervalUtils.getPaddedInterval(startPoint, LINK_DENSITY_WIDTH, dictionary);
            final SVInterval endInterval = SVIntervalUtils.getPaddedInterval(endPoint, LINK_DENSITY_WIDTH, dictionary);
            final int size = startInterval.getLength() + endInterval.getLength() - startInterval.overlapLen(endInterval);
            final Stream<EvidenceTargetLink> startLinks = Utils.stream(links.overlappers(startInterval)).map(SVIntervalTree.Entry::getValue)
                    .filter(link -> link.getPairedStrandedIntervals().getLeft().getInterval().overlaps(startInterval) || link.getPairedStrandedIntervals().getRight().getInterval().overlaps(startInterval));
            final Stream<EvidenceTargetLink> endLinks = Utils.stream(links.overlappers(endInterval)).map(SVIntervalTree.Entry::getValue)
                    .filter(link -> link.getPairedStrandedIntervals().getLeft().getInterval().overlaps(endInterval) || link.getPairedStrandedIntervals().getRight().getInterval().overlaps(endInterval));
            event.linkDensity = Stream.concat(startLinks, endLinks)
                    .distinct()
                    .filter(link -> !eventLinks.contains(link))
                    .mapToInt(link -> link.getReadPairs())
                    .sum() * 1e6 / (double) (size * COPY_NEUTRAL_DEPTH);
        }
    }

    private Map<SimpleSVType.TYPES,List<ReadDepthCluster>> clusterEvents(final SVIntervalTree<LargeSimpleSV> eventsTree, final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector, final SAMSequenceDictionary dictionary) {
        final Set<LargeSimpleSV> visited = new HashSet<>(SVUtils.hashMapCapacity(eventsTree.size()));
        final Map<SimpleSVType.TYPES,List<ReadDepthCluster>> clusteredEvents = new HashMap<>(SVUtils.hashMapCapacity(eventsTree.size()));
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
                    cluster.addAll(Utils.stream(eventsTree.overlappers(clusterEvent.getInterval()))
                            .map(SVIntervalTree.Entry::getValue)
                            .filter(overlapper -> overlapper.getType() == event.getType()) //Cluster by event type
                            .collect(Collectors.toList()));
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
            clusteredEvents.putIfAbsent(event.getType(), new ArrayList<>());
            clusteredEvents.get(event.getType()).add(new ReadDepthCluster(newCluster, copyRatioSegmentOverlapDetector, dictionary));
        }
        return clusteredEvents;
    }

    public Collection<ReadDepthEvent> getEvents() {
        return clusteredEvents.values().stream().flatMap(List::stream).flatMap(entry -> entry.getEventsList().stream()).collect(Collectors.toList());
    }

    private static double metropolisHastingsTemperatureSchedule(final int step) {
        return 1.0;
    }

    public Tuple2<Double,List<ReadDepthEvent>> solve(final JavaSparkContext ctx) {
        List<ReadDepthCluster> clusters = clusteredEvents.values().stream().flatMap(List::stream).collect(Collectors.toList());
        double energy = 0;
        for (int i = 0; i < 1; i++) {
            final JavaRDD<ReadDepthCluster> clusterRdd = ctx.parallelize(clusters);
            final List<Tuple2<Double,ReadDepthCluster>> result = solveStateMaximumPosterior(clusterRdd, parameters, seed + i * 9248837L);
            clusters = result.stream().map(Tuple2::_2).collect(Collectors.toList());
            energy = result.stream().mapToDouble(Tuple2::_1).sum();

            final List<ReadDepthEvent> events = clusters.stream().flatMap(cluster -> cluster.getEventsList().stream()).collect(Collectors.toList());
            final double fractionCalled = events.stream().filter(event -> event.getState() > 0).count() / (double) events.size();

            //solveParameterMaximumPosterior(clusters, parameters, seed + i * 738232L, logger);
            logger.info("Outer iteration " + i + ": f = " + energy + " fraction called = " + fractionCalled);
        }
        final List<ReadDepthEvent> events = clusters.stream().map(ReadDepthCluster::getEventsList).flatMap(List::stream).collect(Collectors.toList());
        for (final ReadDepthCluster cluster : clusters) {
            final double[] eventStates = cluster.getStates();
            for (final ReadDepthEvent event : cluster.getEventsList()) {
                event.copyNumberCallOverlapLikelihood = computeCallOverlapLikelihood(eventStates, event, parameters);
                event.distanceLikelihood = computeCallDistanceLikelihood(eventStates, event, parameters);
                event.copyNumberLikelihood = computeCopyNumberLikelihood(eventStates, event, parameters);
                event.readPairEvidenceLikelihood = computeReadPairEvidenceLikelihood(eventStates, event, parameters);
                event.splitReadEvidenceLikelihood = computeSplitReadEvidenceLikelihood(eventStates, event, parameters);
            }
        }
        return new Tuple2<>(energy,events);
    }

    public Tuple2<Double,List<ReadDepthEvent>> train(final JavaSparkContext ctx, final SVIntervalTree<VariantContext> truthSetTree) {
        List<ReadDepthCluster> clusters = clusteredEvents.values().stream().flatMap(List::stream).collect(Collectors.toList());
        double energy = 0;
        for (int i = 0; i < 1; i++) {
            final JavaRDD<ReadDepthCluster> clusterRdd = ctx.parallelize(clusters);
            final List<Tuple2<Double,ReadDepthCluster>> result = solveStateMaximumPosterior(clusterRdd, parameters, seed + i * 9248837L);
            clusters = result.stream().map(Tuple2::_2).collect(Collectors.toList());
            energy = result.stream().mapToDouble(Tuple2::_1).sum();

            final List<ReadDepthEvent> events = clusters.stream().flatMap(cluster -> cluster.getEventsList().stream()).collect(Collectors.toList());
            final double fractionCalled = events.stream().filter(event -> event.getState() > 0).count() / (double) events.size();

            energy = trainParameters(clusters, parameters, truthSetTree, + i * 738232L, logger);
            logger.info("Training outer iteration " + i + ": f = " + energy + " fraction called = " + fractionCalled);
        }
        final List<ReadDepthEvent> events = clusters.stream().map(ReadDepthCluster::getEventsList).flatMap(List::stream).collect(Collectors.toList());
        for (final ReadDepthCluster cluster : clusters) {
            final double[] eventStates = cluster.getStates();
            for (final ReadDepthEvent event : cluster.getEventsList()) {
                event.copyNumberCallOverlapLikelihood = computeCallOverlapLikelihood(eventStates, event, parameters);
                event.distanceLikelihood = computeCallDistanceLikelihood(eventStates, event, parameters);
                event.copyNumberLikelihood = computeCopyNumberLikelihood(eventStates, event, parameters);
                event.readPairEvidenceLikelihood = computeReadPairEvidenceLikelihood(eventStates, event, parameters);
                event.splitReadEvidenceLikelihood = computeSplitReadEvidenceLikelihood(eventStates, event, parameters);
            }
        }
        setEventTruePositiveFlags(clusters, truthSetTree);
        return new Tuple2<>(energy,events);
    }

    public static double trainParameters(final List<ReadDepthCluster> clusters, final ReadDepthModelParameters parameters, final SVIntervalTree<VariantContext> truthSetTree, final long seed, final Logger logger) {
        final int size = parameters.getParameters().length;
        final NormalDistribution standardNormal = new NormalDistribution(0, 0.1);
        standardNormal.reseedRandomGenerator(seed);
        final Function<double[],Double> energyFunction = x -> trainingEnergyFunction(x, clusters, truthSetTree, seed);
        final Supplier<double[]> sampler = () -> parameterStepSampler(standardNormal);
        final double[] lowerBound = ReadDepthModelParameters.getDefaultLowerBounds();
        final double[] upperBound = ReadDepthModelParameters.getDefaultUpperBounds();
        final double[] x0 = parameters.getParameters();

        //final SimulatedAnnealingSolver parameterSolver = new SimulatedAnnealingSolver(size, energyFunction, sampler, lowerBound, upperBound);
        //final double finalEnergy = parameterSolver.solve(x0, 2000, 100);

        final GradientDescentSolver parameterSolver = new GradientDescentSolver(energyFunction, size, 1e-3, 1e-3);
        final double finalEnergy = parameterSolver.solve(x0, 2000, 100);

        parameters.setParameters(parameterSolver.getSolution());
        logger.info(parameters.toString());
        return finalEnergy;
    }

    public static double trainingEnergyFunction(final double[] x, final List<ReadDepthCluster> clusters, final SVIntervalTree<VariantContext> truthSetTree, final long seed) {
        final ReadDepthModelParameters parametersX = new ReadDepthModelParameters();
        parametersX.setParameters(x);
        final long seedX = seed + new Double(Arrays.stream(x).sum()).hashCode();
        clusters.stream().forEach(cluster -> solveStateMaximumPosterior(cluster, parametersX, seed + seedX));
        return computeTrainingError(clusters, truthSetTree);
    }

    public static void setEventTruePositiveFlags(final List<ReadDepthCluster> eventClusters, final SVIntervalTree<VariantContext> truthSetTree) {
        final List<ReadDepthEvent> calledEvents = eventClusters.stream().flatMap(cluster -> cluster.getEventsList().stream()).filter(event -> event.getState() > 0).collect(Collectors.toList());
        for (final ReadDepthEvent event : calledEvents) {
            event.isTruePositive = SVIntervalUtils.getIntervalsWithReciprocalOverlapInTree(event.getEvent().getInterval(), truthSetTree, 0.5).stream()
                    .anyMatch(entry -> (entry.getValue().getStructuralVariantType() == StructuralVariantType.DEL && event.getEvent().getType() == SimpleSVType.TYPES.DEL)
                            || (entry.getValue().getStructuralVariantType() == StructuralVariantType.DUP && event.getEvent().getType() == SimpleSVType.TYPES.DUP_TAND));
        }
    }

    public static double computeTrainingError(final List<ReadDepthCluster> eventClusters, final SVIntervalTree<VariantContext> truthSetTree) {
        final List<ReadDepthEvent> calledEvents = eventClusters.stream().flatMap(cluster -> cluster.getEventsList().stream()).filter(event -> event.getState() > 0).collect(Collectors.toList());
        final int truePositives = (int) calledEvents.stream()
                .filter(event -> SVIntervalUtils.getIntervalsWithReciprocalOverlapInTree(event.getEvent().getInterval(), truthSetTree, 0.5).stream()
                        .anyMatch(entry -> (entry.getValue().getStructuralVariantType() == StructuralVariantType.DEL && event.getEvent().getType() == SimpleSVType.TYPES.DEL)
                                || (entry.getValue().getStructuralVariantType() == StructuralVariantType.DUP && event.getEvent().getType() == SimpleSVType.TYPES.DUP_TAND)))
                .count();
        final int falsePositives = calledEvents.size() - truePositives;
        return falsePositives - truePositives;
    }

    public static double solveParameterMaximumPosterior(final List<ReadDepthCluster> clusters, final ReadDepthModelParameters parameters, final long seed, final Logger logger) {
        final int size = parameters.getParameters().length;
        final NormalDistribution standardNormal = new NormalDistribution(0, 0.01);
        standardNormal.reseedRandomGenerator(seed);
        final Function<double[],Double> energyFunction = x -> clusters.stream().mapToDouble(cluster -> {
            final ReadDepthModelParameters parametersX = new ReadDepthModelParameters();
            parametersX.setParameters(x);
            final double[] eventStates = cluster.getStates();
            return -computeLogPosterior(eventStates, cluster, parametersX);
        }).sum();
        final Supplier<double[]> sampler = () -> parameterStepSampler(standardNormal);
        final double[] lowerBound = ReadDepthModelParameters.getDefaultLowerBounds();
        final double[] upperBound = ReadDepthModelParameters.getDefaultUpperBounds();
        final SimulatedAnnealingSolver parameterSolver = new SimulatedAnnealingSolver(size, energyFunction, sampler, lowerBound, upperBound);

        final double[] x0 = parameters.getParameters();
        final double finalEnergy = parameterSolver.solve(x0, 10000, 0);
        parameters.setParameters(parameterSolver.getSolution());
        logger.info(parameters.toString());
        return finalEnergy;
    }

    private static List<Tuple2<Double,ReadDepthCluster>> solveStateMaximumPosterior(final JavaRDD<ReadDepthCluster> clusters, final ReadDepthModelParameters parameters, final long seed) {
        return clusters.map(cluster -> solveStateMaximumPosterior(cluster, parameters, seed)).collect();
    }

    private static double[] stateStepSampler(final int size, final Random random) {
        /*final double[] sample = standardNormal.sample(size);
        for (int i = 0; i < sample.length; i++) {
            sample[i] *= 0.1;
        }
        return sample;*/
        final double[] sample = new double[size];
        boolean isZero = true;
        while (isZero) {
            for (int i = 0; i < size; i++) {
                sample[i] = random.nextInt(3) - 1;
                if (sample[i] != 0) {
                    isZero = false;
                }
            }
        }
        return sample;
    }

    public static Tuple2<Double,ReadDepthCluster> solveStateMaximumPosterior(final ReadDepthCluster cluster, final ReadDepthModelParameters parameters, final long seed) {
        if (cluster.getEventsList().size() == 1) {
            return solveStateMaximumPosteriorBF(cluster, parameters);
        }
        return solveStateMaximumPosteriorSA(cluster, parameters, seed);
    }

    private static Tuple2<Double,ReadDepthCluster> solveStateMaximumPosteriorBF(final ReadDepthCluster cluster, final ReadDepthModelParameters parameters) {
        final int numEvents = cluster.getEventsTree().size();
        if (numEvents != 1) {
            throw new GATKException("Brute force only supports one event");
        }
        final double[] x = new double[numEvents];
        double minX = 0;
        double minEnergy = -computeLogPosterior(x, cluster, parameters);
        final int maxState = (int) (cluster.getType() == DEL ? MAX_PLOIDY : MAX_PLOIDY * MAX_DUPLICATION_NUMBER);
        for (int i = 1; i <= maxState; i++) {
            x[0] = i;
            final double energy = -computeLogPosterior(x, cluster, parameters);
            if (i == 0 || energy < minEnergy) {
                minX = i;
                minEnergy = energy;
            }
        }
        cluster.getEventsList().get(0).setState(minX);
        return new Tuple2<>(minEnergy, cluster);
    }

    private static Tuple2<Double,ReadDepthCluster> solveStateMaximumPosteriorSA(final ReadDepthCluster cluster, final ReadDepthModelParameters parameters, final long seed) {
        final List<ReadDepthEvent> events = cluster.getEventsList();
        final Function<double[], Double> energyFunction = x -> -computeLogPosterior(x, cluster, parameters);
        final int size = events.size();
        final Random random = new Random(seed);
        final Supplier<double[]> sampler = () -> stateStepSampler(size, random);
        final double[] lowerBound = new double[size];
        final double[] upperBound = new double[size];
        Arrays.fill(lowerBound, 0);
        final int maxState = (int) (cluster.getType() == DEL ? MAX_PLOIDY : MAX_PLOIDY * MAX_DUPLICATION_NUMBER);
        Arrays.fill(upperBound, maxState);
        final SimulatedAnnealingSolver sa = new SimulatedAnnealingSolver(size, energyFunction, sampler, lowerBound, upperBound);

        final double[] x0 = new double[size];
        for (int i = 0; i < events.size(); i++) {
            x0[i] = events.get(i).getState();
        }
        final int numSteps = (int) Math.min(10000, Math.pow(10, size));
        final double finalEnergy = sa.solve(x0, numSteps, 0);
        //final double finalEnergy = sa.solve(x0, numSteps, ReadDepthModel::metropolisHastingsTemperatureSchedule, 100);

        final double[] minEnergyState = sa.getSolution();
        for (int i = 0; i < events.size(); i++) {
            events.get(i).setState(minEnergyState[i]);
        }
        return new Tuple2<>(finalEnergy, cluster);
    }

    private static double computeLogPosterior(final double[] x, final ReadDepthCluster cluster, final ReadDepthModelParameters parameters) {
        return computeLogLikelihood(x, cluster, parameters) + computeLogPrior(x, cluster, parameters);
    }

    private static double computeLogLikelihood(final double[] x, final ReadDepthCluster cluster, final ReadDepthModelParameters parameters) {
        return cluster.getEventsList().stream().mapToDouble(event -> computeCopyNumberLikelihood(x, event, parameters)
                + computeReadPairEvidenceLikelihood(x, event, parameters)
                + computeSplitReadEvidenceLikelihood(x, event, parameters)
                + computeCallDistanceLikelihood(x, event, parameters)
                //+ computeLinkDensityLikelihood(x, event, parameters)
                + computeCallOverlapLikelihood(x, event, parameters)).sum();
    }

    private static double computeLogPrior(final double[] x, final ReadDepthCluster cluster, final ReadDepthModelParameters parameters) {
        return cluster.getEventsList().stream().mapToDouble(event -> computePloidyPrior(x, event, cluster.getType(), parameters)).sum();
    }

    private static double computeLinkDensityLikelihood(final double[] x, final ReadDepthEvent event, final ReadDepthModelParameters parameters) {
        final int id = event.getId();
        final double mean = x[id] == 0 ? 30 : 5;
        final double std = x[id] == 0 ? 10 : 5;
        return unscaledLogNormal(event.linkDensity, mean, std);
    }

    private static double computeCallOverlapLikelihood(final double[] x, final ReadDepthEvent event, final ReadDepthModelParameters parameters) {
        final int id = event.getId();
        final double mean = x[id] == 0 ? parameters.getParameter(ReadDepthModelParameters.ParameterEnum.CALL_OVERLAP_MEAN_0) : parameters.getParameter(ReadDepthModelParameters.ParameterEnum.CALL_OVERLAP_MEAN_1);
        final double std = x[id] == 0 ? parameters.getParameter(ReadDepthModelParameters.ParameterEnum.CALL_OVERLAP_STD_0) : parameters.getParameter(ReadDepthModelParameters.ParameterEnum.CALL_OVERLAP_STD_1);
        return unscaledLogNormal(event.copyNumberCallOverlap, mean, std);
    }

    private static double computeCallDistanceLikelihood(final double[] x, final ReadDepthEvent event, final ReadDepthModelParameters parameters) {
        final int id = event.getId();
        final double mean = x[id] == 0 ? parameters.getParameter(ReadDepthModelParameters.ParameterEnum.CALL_START_DISTANCE_MEAN_0) : parameters.getParameter(ReadDepthModelParameters.ParameterEnum.CALL_START_DISTANCE_MEAN_1);
        final double std = x[id] == 0 ? parameters.getParameter(ReadDepthModelParameters.ParameterEnum.CALL_START_DISTANCE_STD_0) : parameters.getParameter(ReadDepthModelParameters.ParameterEnum.CALL_START_DISTANCE_STD_1);
        return unscaledLogNormal(event.leftDistance, mean, std) + unscaledLogNormal(event.rightDistance, mean, std);
    }

    private static double computeReadPairEvidenceLikelihood(final double[] x, final ReadDepthEvent event, final ReadDepthModelParameters parameters) {
        final double singleCopyDepth = COPY_NEUTRAL_DEPTH * 0.5;
        final double expectedReadEvidenceFraction = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.EXPECTED_READ_PAIR_EVIDENCE_FRACTION);
        final double expectedReadEvidenceStd = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.EXPECTED_READ_PAIR_EVIDENCE_STD);
        final double expectedEvidence = x[event.getId()] * singleCopyDepth * expectedReadEvidenceFraction;
        final double sigma = expectedReadEvidenceStd * singleCopyDepth;
        return unscaledLogNormal(event.getEvent().getReadPairEvidence(), expectedEvidence, sigma);
    }

    private static double computeSplitReadEvidenceLikelihood(final double[] x, final ReadDepthEvent event, final ReadDepthModelParameters parameters) {
        final double singleCopyDepth = COPY_NEUTRAL_DEPTH * 0.5;
        final double expectedReadEvidenceFraction = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.EXPECTED_SPLIT_READ_EVIDENCE_FRACTION);
        final double expectedReadEvidenceStd = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.EXPECTED_SPLIT_READ_EVIDENCE_STD);
        final double expectedEvidence = x[event.getId()] * singleCopyDepth * expectedReadEvidenceFraction;
        final double sigma = expectedReadEvidenceStd * singleCopyDepth;
        return unscaledLogNormal(event.getEvent().getSplitReadEvidence(), expectedEvidence, sigma);
    }

    private static double computeCopyNumberLikelihood(final double[] x, final ReadDepthEvent event, final ReadDepthModelParameters parameters) {
        final double dataCopyNumber = event.observedCopyNumber + event.optimizedOverlapInfoList.stream().mapToDouble(pair -> x[pair._1] * pair._2).sum();
        /*final double expectedCopyNumber = event.overlapInfoList.stream()
                .mapToDouble(overlapInfo -> overlapInfo.weight * (overlapInfo.segmentCopyNumber - overlapInfo.overlappingIds.stream()
                        .mapToDouble(id -> x[id])
                        .sum()))
                .sum();*/
        final double modelCopyNumber = 2 + x[event.getId()] * (event.getEvent().getType() == DEL ? -1 : 1);
        final double sigma = modelCopyNumber == 0 ? parameters.getParameter(ReadDepthModelParameters.ParameterEnum.COPY_NUMBER_STD_0) : parameters.getParameter(ReadDepthModelParameters.ParameterEnum.COPY_NUMBER_STD_1);
        return unscaledLogNormal(dataCopyNumber, modelCopyNumber, sigma);
    }

    private static double computePloidyPrior(final double[] x, final ReadDepthEvent event,
                                      final SimpleSVType.TYPES type, final ReadDepthModelParameters parameters) {
        final List<Integer> overlapperIds = event.overlappingEventIds;
        final double parameterConstraintStd = parameters.getParameter(ReadDepthModelParameters.ParameterEnum.PARAMETER_CONSTRAINT_STD);
        double totalPloidy = overlapperIds.stream().mapToDouble(id -> getPloidy(x[id], type)).sum();
        final double excessPloidy;
        if (totalPloidy > MAX_PLOIDY) {
            excessPloidy = totalPloidy - MAX_PLOIDY;
        } else {
            excessPloidy = 0;
        }
        return unscaledLogNormal(excessPloidy, 0, parameterConstraintStd);
    }

    private static double getPloidy(final double val, final SimpleSVType.TYPES type) {
        return type == SimpleSVType.TYPES.DUP_TAND ? Math.min(val, 1) : val;
    }

    private final static class ReadDepthModelParameters implements Serializable {

        public static final long serialVersionUID = 1L;
        public enum ParameterEnum {
            PARAMETER_CONSTRAINT_STD,
            EXPECTED_READ_PAIR_EVIDENCE_FRACTION,
            EXPECTED_READ_PAIR_EVIDENCE_STD,
            EXPECTED_SPLIT_READ_EVIDENCE_FRACTION,
            EXPECTED_SPLIT_READ_EVIDENCE_STD,
            COPY_NUMBER_STD_0,
            COPY_NUMBER_STD_1,
            CALL_START_DISTANCE_MEAN_0,
            CALL_START_DISTANCE_STD_0,
            CALL_START_DISTANCE_MEAN_1,
            CALL_START_DISTANCE_STD_1,
            CALL_OVERLAP_MEAN_0,
            CALL_OVERLAP_STD_0,
            CALL_OVERLAP_MEAN_1,
            CALL_OVERLAP_STD_1,
        }

        public static final double DEFAULT_PARAMETER_CONSTRAINT_STD = 0.001;
        public static final double DEFAULT_EXPECTED_READ_PAIR_EVIDENCE_FRACTION = 0.5;
        public static final double DEFAULT_EXPECTED_READ_PAIR_EVIDENCE_STD = 0.5;
        public static final double DEFAULT_EXPECTED_SPLIT_READ_EVIDENCE_FRACTION = 0.5;
        public static final double DEFAULT_EXPECTED_SPLIT_READ_EVIDENCE_STD = 1;
        public static final double DEFAULT_COPY_NUMBER_STD_0 = 1;
        public static final double DEFAULT_COPY_NUMBER_STD_1 = 1;
        public static final double DEFAULT_CALL_START_DISTANCE_MEAN_0 = 6.5;
        public static final double DEFAULT_CALL_START_DISTANCE_STD_0 = 1.3;
        public static final double DEFAULT_CALL_START_DISTANCE_MEAN_1 = 4.2;
        public static final double DEFAULT_CALL_START_DISTANCE_STD_1 = 1.4;
        public static final double DEFAULT_CALL_OVERLAP_MEAN_0 = 0.5;
        public static final double DEFAULT_CALL_OVERLAP_STD_0 = 0.1;
        public static final double DEFAULT_CALL_OVERLAP_MEAN_1 = 1;
        public static final double DEFAULT_CALL_OVERLAP_STD_1 = 0.05;

        private double[] parameters;

        public ReadDepthModelParameters() {
            parameters = new double[ParameterEnum.values().length];
            setParameter(ParameterEnum.PARAMETER_CONSTRAINT_STD, DEFAULT_PARAMETER_CONSTRAINT_STD);
            setParameter(ParameterEnum.EXPECTED_READ_PAIR_EVIDENCE_FRACTION, DEFAULT_EXPECTED_READ_PAIR_EVIDENCE_FRACTION);
            setParameter(ParameterEnum.EXPECTED_READ_PAIR_EVIDENCE_STD, DEFAULT_EXPECTED_READ_PAIR_EVIDENCE_STD);
            setParameter(ParameterEnum.EXPECTED_SPLIT_READ_EVIDENCE_FRACTION, DEFAULT_EXPECTED_SPLIT_READ_EVIDENCE_FRACTION);
            setParameter(ParameterEnum.EXPECTED_SPLIT_READ_EVIDENCE_STD, DEFAULT_EXPECTED_SPLIT_READ_EVIDENCE_STD);
            setParameter(ParameterEnum.COPY_NUMBER_STD_0, DEFAULT_COPY_NUMBER_STD_0);
            setParameter(ParameterEnum.COPY_NUMBER_STD_1, DEFAULT_COPY_NUMBER_STD_1);
            setParameter(ParameterEnum.CALL_START_DISTANCE_MEAN_0, DEFAULT_CALL_START_DISTANCE_MEAN_0);
            setParameter(ParameterEnum.CALL_START_DISTANCE_STD_0, DEFAULT_CALL_START_DISTANCE_STD_0);
            setParameter(ParameterEnum.CALL_START_DISTANCE_MEAN_1, DEFAULT_CALL_START_DISTANCE_MEAN_1);
            setParameter(ParameterEnum.CALL_START_DISTANCE_STD_1, DEFAULT_CALL_START_DISTANCE_STD_1);
            setParameter(ParameterEnum.CALL_OVERLAP_MEAN_0, DEFAULT_CALL_OVERLAP_MEAN_0);
            setParameter(ParameterEnum.CALL_OVERLAP_STD_0, DEFAULT_CALL_OVERLAP_STD_0);
            setParameter(ParameterEnum.CALL_OVERLAP_MEAN_1, DEFAULT_CALL_OVERLAP_MEAN_1);
            setParameter(ParameterEnum.CALL_OVERLAP_STD_1, DEFAULT_CALL_OVERLAP_STD_1);
        }

        public double[] getParameters() {
            return parameters;
        }

        public void setParameter(final ParameterEnum parameter, final double value) {
            parameters[parameter.ordinal()] = value;
        }

        public double getParameter(final ParameterEnum parameter) {
            return parameters[parameter.ordinal()];
        }

        public void setParameters(final double[] values) {
            Utils.nonNull(values, "Parameter vector cannot be null");
            Utils.validateArg(values.length == ParameterEnum.values().length, "Invalid parameter vector size");
            parameters = values;
        }

        public static double[] getDefaultUpperBounds() {
            final double[] params = new double[ParameterEnum.values().length];
            params[ParameterEnum.PARAMETER_CONSTRAINT_STD.ordinal()] = 0.1;
            params[ParameterEnum.EXPECTED_READ_PAIR_EVIDENCE_FRACTION.ordinal()] = 2;
            params[ParameterEnum.EXPECTED_READ_PAIR_EVIDENCE_STD.ordinal()] = 10;
            params[ParameterEnum.EXPECTED_SPLIT_READ_EVIDENCE_FRACTION.ordinal()] = 2;
            params[ParameterEnum.EXPECTED_SPLIT_READ_EVIDENCE_STD.ordinal()] = 10;
            params[ParameterEnum.COPY_NUMBER_STD_0.ordinal()] = 10;
            params[ParameterEnum.COPY_NUMBER_STD_1.ordinal()] = 10;
            params[ParameterEnum.CALL_START_DISTANCE_MEAN_0.ordinal()] = 100;
            params[ParameterEnum.CALL_START_DISTANCE_STD_0.ordinal()] = 100;
            params[ParameterEnum.CALL_START_DISTANCE_MEAN_1.ordinal()] = 100;
            params[ParameterEnum.CALL_START_DISTANCE_STD_1.ordinal()] = 100;
            params[ParameterEnum.CALL_OVERLAP_MEAN_0.ordinal()] = 1;
            params[ParameterEnum.CALL_OVERLAP_STD_0.ordinal()] = 100;
            params[ParameterEnum.CALL_OVERLAP_MEAN_1.ordinal()] = 1;
            params[ParameterEnum.CALL_OVERLAP_STD_1.ordinal()] = 100;
            return params;
        }

        public static double[] getDefaultLowerBounds() {
            final double[] params = new double[ParameterEnum.values().length];
            params[ParameterEnum.PARAMETER_CONSTRAINT_STD.ordinal()] = 1e-3;
            params[ParameterEnum.EXPECTED_READ_PAIR_EVIDENCE_FRACTION.ordinal()] = 0;
            params[ParameterEnum.EXPECTED_READ_PAIR_EVIDENCE_STD.ordinal()] = 1e-3;
            params[ParameterEnum.EXPECTED_SPLIT_READ_EVIDENCE_FRACTION.ordinal()] = 0;
            params[ParameterEnum.EXPECTED_SPLIT_READ_EVIDENCE_STD.ordinal()] = 1e-3;
            params[ParameterEnum.COPY_NUMBER_STD_0.ordinal()] = 1e-3;
            params[ParameterEnum.COPY_NUMBER_STD_1.ordinal()] = 1e-3;
            params[ParameterEnum.CALL_START_DISTANCE_MEAN_0.ordinal()] = 0;
            params[ParameterEnum.CALL_START_DISTANCE_STD_0.ordinal()] = 1e-3;
            params[ParameterEnum.CALL_START_DISTANCE_MEAN_1.ordinal()] = 0;
            params[ParameterEnum.CALL_START_DISTANCE_STD_1.ordinal()] = 1e-3;
            params[ParameterEnum.CALL_OVERLAP_MEAN_0.ordinal()] = 0;
            params[ParameterEnum.CALL_OVERLAP_STD_0.ordinal()] = 1e-3;
            params[ParameterEnum.CALL_OVERLAP_MEAN_1.ordinal()] = 0;
            params[ParameterEnum.CALL_OVERLAP_STD_1.ordinal()] = 1e-3;
            return params;
        }

        @Override
        public String toString() {
            final StringBuilder stringBuilder = new StringBuilder();
            for (final ParameterEnum type : ParameterEnum.values()) {
                stringBuilder.append("{");
                stringBuilder.append(type.toString());
                stringBuilder.append(": ");
                stringBuilder.append(getParameter(type));
                stringBuilder.append("}");
            }
            return stringBuilder.toString();
        }
    }

    final static double[] parameterStepSampler(final NormalDistribution standardNormal) {
        final double[] sample = standardNormal.sample(ReadDepthModelParameters.ParameterEnum.values().length);
        sample[ReadDepthModelParameters.ParameterEnum.PARAMETER_CONSTRAINT_STD.ordinal()] = 0;
        sample[ReadDepthModelParameters.ParameterEnum.EXPECTED_READ_PAIR_EVIDENCE_FRACTION.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.EXPECTED_READ_PAIR_EVIDENCE_STD.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.EXPECTED_SPLIT_READ_EVIDENCE_FRACTION.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.EXPECTED_SPLIT_READ_EVIDENCE_STD.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.COPY_NUMBER_STD_0.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.COPY_NUMBER_STD_1.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.CALL_START_DISTANCE_MEAN_0.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.CALL_START_DISTANCE_STD_0.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.CALL_START_DISTANCE_MEAN_1.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.CALL_START_DISTANCE_STD_1.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.CALL_OVERLAP_MEAN_0.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.CALL_OVERLAP_STD_0.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.CALL_OVERLAP_MEAN_1.ordinal()] *= 1;
        sample[ReadDepthModelParameters.ParameterEnum.CALL_OVERLAP_STD_1.ordinal()] *= 1;
        return sample;
    }

    static final class OverlapInfo implements Serializable {
        public static final long serialVersionUID = 1L;
        public List<Tuple2<Integer,Integer>> overlappingIdsAndSigns;
        public double weight;
        public double segmentCopyNumber;

        public OverlapInfo(List<Tuple2<Integer,Integer>> overlappingIds, double weight, double segmentCopyNumber) {
            this.overlappingIdsAndSigns = overlappingIds;
            this.weight = weight;
            this.segmentCopyNumber = segmentCopyNumber;
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
