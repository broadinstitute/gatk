package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import org.broadinstitute.hellbender.tools.spark.sv.DiscoverVariantsFromReadDepthArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Uses an SV graph and copy number posteriors to call SVs
 */
public final class ReadDepthSVCaller {

    private final SVGraph graph;
    private final DiscoverVariantsFromReadDepthArgumentCollection arguments;
    private final SVIntervalTree<SVCopyNumberInterval> copyNumberPosteriorsTree;

    public ReadDepthSVCaller(final SVGraph graph, final SVIntervalTree<SVCopyNumberInterval> copyNumberPosteriorsTree, final DiscoverVariantsFromReadDepthArgumentCollection arguments) {
        this.graph = graph;
        this.arguments = arguments;
        this.copyNumberPosteriorsTree = copyNumberPosteriorsTree;
    }

    /**
     * Returns true if the two edges each partially overlap one another or if they meet minimum reciprocal overlap
     */
    private static boolean graphPartitioningFunction(final SVInterval a, final SVInterval b, final double minReciprocalOverlap) {
        if (!a.overlaps(b)) return false;
        return (a.getStart() <= b.getStart() && a.getEnd() <= b.getEnd())
                || (b.getStart() <= a.getStart() && b.getEnd() <= a.getEnd())
                || SVIntervalUtils.hasReciprocalOverlap(a, b, minReciprocalOverlap);
    }

    /**
     * Main function that enumerates graph paths, integrates event probabilities, and produces calls
     */
    public static Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> generateEvents(final SVGraph graph, final int groupId, final double minEventProb,
                                                                                                           final double maxPathLengthFactor, final int maxEdgeVisits,
                                                                                                           final SVIntervalTree<SVCopyNumberInterval> copyNumberPosteriorsTree,
                                                                                                           final int maxQueueSize, final int ploidy, final int minSize, final double minHaplotypeProb) {
        final SVGraphGenotyper searcher = new SVGraphGenotyper(graph);
        final Collection<IndexedSVGraphPath> paths = searcher.enumerate(maxPathLengthFactor, maxEdgeVisits, maxQueueSize);
        if (paths == null) return null;
        final Collection<SVGraphGenotype> genotypes = enumerateGenotypes(paths, graph, copyNumberPosteriorsTree, groupId);
        setDepthProbabilities(genotypes);
        setGenotypeEvidenceProbabilities(genotypes, graph);
        setProbabilities(genotypes);
        final Collection<Tuple2<CalledSVGraphEvent, Double>> integratedEvents = integrateEdgeEvents(genotypes, graph);
        final Collection<CalledSVGraphEvent> probabilityFilteredEvents = filterEventsByProbability(integratedEvents, minEventProb);
        final Collection<CalledSVGraphEvent> mergedEvents = probabilityFilteredEvents; //TODO mergeAdjacentEvents(probabilityFilteredEvents);
        final Collection<CalledSVGraphEvent> sizeFilteredEvents = filterEventsBySize(mergedEvents, minSize);
        final Collection<CalledSVGraphGenotype> haplotypes = convertToCalledHaplotypes(filterHaplotypesByProbability(genotypes, minHaplotypeProb), graph);
        return new Tuple2<>(haplotypes, sizeFilteredEvents);
    }

    private static Collection<CalledSVGraphGenotype> convertToCalledHaplotypes(final Collection<SVGraphGenotype> haplotypes, final SVGraph graph) {
        return haplotypes.stream().map(h -> new CalledSVGraphGenotype(h, graph)).collect(Collectors.toList());
    }

    private static Collection<SVGraphGenotype> filterHaplotypesByProbability(final Collection<SVGraphGenotype> haplotypes, final double minProb) {
        return haplotypes.stream().filter(h -> h.getDepthProbability() >= minProb).collect(Collectors.toList());
    }

    private static Collection<CalledSVGraphEvent> filterEventsByProbability(final Collection<Tuple2<CalledSVGraphEvent, Double>> calledSVGraphEvents, final double minProb) {
        return calledSVGraphEvents.stream().filter(pair -> pair._2 >= minProb).map(Tuple2::_1).collect(Collectors.toList());
    }

    private static Collection<CalledSVGraphEvent> filterEventsBySize(final Collection<CalledSVGraphEvent> calledSVGraphEvents, final int minSize) {
        return calledSVGraphEvents.stream().filter(sv -> sv.getInterval().getLength() >= minSize).collect(Collectors.toList());
    }

    private static Collection<SVGraphGenotype> enumerateGenotypes(final Collection<IndexedSVGraphPath> paths,
                                                                  final SVGraph graph,
                                                                  final SVIntervalTree<SVCopyNumberInterval> copyNumberPosteriorsTree,
                                                                  final int groupId) {
        final int numEdges = graph.getEdges().size();
        final List<EdgeCopyNumberPosterior> copyNumberPosteriors = getEdgeCopyNumberPosteriors(copyNumberPosteriorsTree, graph);
        //TODO: assumes diploid
        int genotypeId = 0;
        if (paths.size() > 30000) return Collections.emptyList(); //Genotypes wouldn't fit into an array
        final List<IndexedSVGraphPath> pathsList = new ArrayList<>(paths);
        final Collection<SVGraphGenotype> genotypes = new ArrayList<>(paths.size() * paths.size());
        final List<int[]> edgeCopyNumberStates = pathsList.stream().map(path -> getEdgeCopyNumberStates(path, numEdges)).collect(Collectors.toList());
        for (int i = 0; i < pathsList.size(); i++) {
            final IndexedSVGraphPath pathA = pathsList.get(i);
            final int[] statesA = edgeCopyNumberStates.get(i);
            for (int j = i; j < pathsList.size(); j++) {
                final IndexedSVGraphPath pathB = pathsList.get(j);
                final int[] statesB = edgeCopyNumberStates.get(j);

                final Collection<IndexedSVGraphPath> genotypePaths = new ArrayList<>(2);
                genotypePaths.add(pathA);
                genotypePaths.add(pathB);

                final int[] genotypeStates = new int[numEdges];
                for (int k = 0; k < numEdges; k++) {
                    genotypeStates[k] = statesA[k] + statesB[k];
                }
                final double logLikelihood = IntStream.range(0, copyNumberPosteriors.size()).mapToDouble(k -> copyNumberPosteriors.get(k).getLogPosterior(genotypeStates[k])).sum();

                final SVGraphGenotype genotype = new SVGraphGenotype(groupId, genotypeId, genotypePaths);
                genotype.setDepthLikelihood(logLikelihood);
                genotypes.add(genotype);
                genotypeId++;
            }
        }
        return genotypes;
    }

    private static int[] getEdgeCopyNumberStates(final IndexedSVGraphPath path, final int numEdges) {
        final int[] states = new int[numEdges];
        for (final IndexedSVGraphEdge edge : path.getEdges()) {
            states[edge.getIndex()]++;
        }
        return states;
    }



    /**
     * Gets copy number posteriors for all edges
     */
    private static List<EdgeCopyNumberPosterior> getEdgeCopyNumberPosteriors(final SVIntervalTree<SVCopyNumberInterval> copyNumberPosteriorsTree,
                                                                             final SVGraph graph) {

        final List<IndexedSVGraphEdge> edges = graph.getEdges();
        final List<IndexedSVGraphEdge> referenceEdges = graph.getReferenceEdges();
        final List<SVCopyNumberInterval> copyNumberIntervals = getSortedOverlappingCopyNumberIntervals(referenceEdges, copyNumberPosteriorsTree);
        if (copyNumberIntervals.isEmpty()) return Collections.emptyList();

        final SVIntervalTree<SVCopyNumberInterval> copyNumberIntervalTree = SVIntervalUtils.buildCopyNumberIntervalTree(copyNumberIntervals);
        final List<EdgeCopyNumberPosterior> edgeCopyNumberPosteriors = new ArrayList<>(edges.size());

        final int numCopyNumberStates = copyNumberIntervals.get(0).getCopyNumberLogPosteriorsArray().length;
        for (final SVCopyNumberInterval copyNumberInterval : copyNumberIntervals) {
            Utils.validate(copyNumberInterval.getCopyNumberLogPosteriorsArray().length == numCopyNumberStates, "Dimension of copy number interval posteriors is not consistent");
        }

        for (final IndexedSVGraphEdge edge : edges) {
            final double[] edgePosterior = computeEdgePosteriors(edge, copyNumberIntervalTree, numCopyNumberStates);
            edgeCopyNumberPosteriors.add(new EdgeCopyNumberPosterior(edgePosterior));
        }

        return edgeCopyNumberPosteriors;
    }

    /**
     * Calculates (approximate) copy number posterior likelihoods for the given edge
     */
    private static final double[] computeEdgePosteriors(final IndexedSVGraphEdge edge, final SVIntervalTree<SVCopyNumberInterval> copyNumberIntervalTree, final int numCopyNumberStates) {
        final double[] copyNumberPosterior = new double[numCopyNumberStates];
        if (edge.isReference()) {
            final SVInterval edgeInterval = edge.getInterval();
            final List<Tuple2<double[], SVInterval>> posteriorsAndIntervalsList = Utils.stream(copyNumberIntervalTree.overlappers(edgeInterval))
                    .map(entry -> new Tuple2<>(entry.getValue().getCopyNumberLogPosteriorsArray(), entry.getInterval()))
                    .collect(Collectors.toList());
            for (final Tuple2<double[], SVInterval> posteriorAndInterval : posteriorsAndIntervalsList) {
                final double[] intervalPosterior = posteriorAndInterval._1;
                final SVInterval interval = posteriorAndInterval._2;
                //Down-weights partially-overlapping intervals
                final double weight = interval.overlapLen(edgeInterval) / (double) interval.getLength();
                for (int i = 0; i < numCopyNumberStates; i++) {
                    copyNumberPosterior[i] += intervalPosterior[i] * weight;
                }
            }
        }
        return copyNumberPosterior;
    }

    private static List<SVCopyNumberInterval> getSortedOverlappingCopyNumberIntervals(final Collection<IndexedSVGraphEdge> edges, final SVIntervalTree<SVCopyNumberInterval> copyNumberPosteriorsTree) {
        return edges.stream()
                .flatMap(edge -> Utils.stream(copyNumberPosteriorsTree.overlappers(edge.getInterval())))
                .map(SVIntervalTree.Entry::getValue)
                .distinct()
                .sorted(SVIntervalUtils.getCopyNumberIntervalDictionaryOrderComparator())
                .collect(Collectors.toList());
    }

    /**
     * Normalizes genotype likelihoods
     */
    private static final void setDepthProbabilities(final Collection<SVGraphGenotype> genotypes) {
        final double logDenom = Math.log(genotypes.stream().mapToDouble(SVGraphGenotype::getDepthLikelihood)
                .map(Math::exp)
                .map(val -> Math.max(val, Double.MIN_NORMAL))
                .sum());
        for (final SVGraphGenotype genotype : genotypes) {
            genotype.setDepthProbability(Math.exp(genotype.getDepthLikelihood() - logDenom));
        }
    }

    private static final void setProbabilities(final Collection<SVGraphGenotype> genotypes) {
        final double denom = genotypes.stream().mapToDouble(g -> g.getEvidenceProbability() * g.getDepthProbability()).sum();
        for (final SVGraphGenotype genotype : genotypes) {
            genotype.setProbability(genotype.getEvidenceProbability() * genotype.getDepthProbability() / denom);
        }
    }

    /**
     * Counts the number of times each reference edge was visited and whether each was inverted
     */
    private static Tuple2<List<int[]>, List<boolean[]>> getReferenceEdgeCountsAndInversions(final SVGraphGenotype haplotypes, final List<IndexedSVGraphEdge> edges) {
        final int numHaplotypes = haplotypes.getHaplotypes().size();
        final List<int[]> referenceEdgeCountsList = new ArrayList<>(numHaplotypes);
        final List<boolean[]> referenceEdgeInversionsList = new ArrayList<>(numHaplotypes);
        for (int haplotypeId = 0; haplotypeId < numHaplotypes; haplotypeId++) {
            final int[] referenceEdgeCounts = new int[edges.size()];
            final boolean[] referenceEdgeInversions = new boolean[edges.size()];
            final IndexedSVGraphPath path = haplotypes.getHaplotypes().get(haplotypeId);
            for (final IndexedSVGraphEdge edge : path.getEdges()) {
                final int edgeIndex = edge.getIndex();
                if (edge.isReference()) {
                    referenceEdgeCounts[edgeIndex]++;
                    if (edge.isInverted()) {
                        referenceEdgeInversions[edgeIndex] = true; //True when visited at least once
                    }
                }
            }
            referenceEdgeCountsList.add(referenceEdgeCounts);
            referenceEdgeInversionsList.add(referenceEdgeInversions);
        }
        return new Tuple2<>(referenceEdgeCountsList, referenceEdgeInversionsList);
    }

    /**
     * Creates events occurring at the given reference edge
     */
    private static final List<SVGraphEvent> getEdgeEvents(final IndexedSVGraphEdge edge,
                                                          final int groupId,
                                                          final int pathId,
                                                          final double probability,
                                                          final double evidenceProbability,
                                                          final List<int[]> referenceEdgeCountsList,
                                                          final List<boolean[]> referenceEdgeInversionsList,
                                                          final List<SVGraphNode> nodes) {
        final int edgeIndex = edge.getIndex();
        boolean deletion = referenceEdgeCountsList.stream().anyMatch(counts -> counts[edgeIndex] < 1);
        boolean duplication = referenceEdgeCountsList.stream().anyMatch(counts -> counts[edgeIndex] > 1);
        boolean inversion = referenceEdgeInversionsList.stream().anyMatch(arr -> arr[edgeIndex]);
        final int numEvents = (deletion ? 1 : 0) + (duplication ? 1 : 0) + (inversion ? 1 : 0) - (deletion && duplication ? 1 : 0);
        final List<SVGraphEvent> events = new ArrayList<>(numEvents);
        if (numEvents > 0) {
            final SVGraphNode nodeA = nodes.get(edge.getNodeAIndex());
            final SVGraphNode nodeB = nodes.get(edge.getNodeBIndex());
            final int start = nodeA.getPosition();
            final int end = nodeB.getPosition();
            final SVInterval interval = new SVInterval(nodeA.getContig(), start, end);
            if (deletion) {
                events.add(new SVGraphEvent(CalledSVGraphEvent.Type.DEL, interval, groupId, pathId, probability, evidenceProbability, true));
            }
            if (duplication && inversion) {
                events.add(new SVGraphEvent(CalledSVGraphEvent.Type.DUP_INV, interval, groupId, pathId, probability, evidenceProbability, true));
            } else {
                if (duplication) {
                    events.add(new SVGraphEvent(CalledSVGraphEvent.Type.DUP, interval, groupId, pathId, probability, evidenceProbability, true));
                }
                if (inversion) {
                    events.add(new SVGraphEvent(CalledSVGraphEvent.Type.INV, interval, groupId, pathId, probability, evidenceProbability, true));
                }
            }
        }
        return events;
    }

    /**
     * Aggregates potential event calls for the given haplotypes, producing a map from reference edge index to list of events
     */
    private static Map<Integer, List<SVGraphEvent>> getHaplotypeEvents(final SVGraphGenotype haplotypes, final SVGraph graph, final IndexedSVGraphPath referencePath) {
        final List<IndexedSVGraphEdge> edges = graph.getEdges();
        final Tuple2<List<int[]>, List<boolean[]>> referenceEdgeResults = getReferenceEdgeCountsAndInversions(haplotypes, edges);
        final List<int[]> referenceEdgeCountsList = referenceEdgeResults._1;
        final List<boolean[]> referenceEdgeInversionsList = referenceEdgeResults._2;

        final List<SVGraphNode> nodes = graph.getNodes();
        final Map<Integer, List<SVGraphEvent>> events = new HashMap<>(SVUtils.hashMapCapacity((edges).size()));
        final int groupId = haplotypes.getGroupId();
        final int pathId = haplotypes.getGenotypeId();
        final double probability = haplotypes.getDepthProbability();
        final double evidenceProbability = haplotypes.getEvidenceProbability();
        for (final IndexedSVGraphEdge edge : referencePath.getEdges()) {
            events.put(edge.getIndex(), getEdgeEvents(edge, groupId, pathId, probability, evidenceProbability, referenceEdgeCountsList, referenceEdgeInversionsList, nodes));
        }
        return events;
    }

    private static double sumProbabilities(final Collection<SVGraphEvent> events) {
        return events.stream().mapToDouble(SVGraphEvent::getProbability).sum();
    }

    private static <T> Map<Integer, List<T>> flattenMaps(final Collection<Map<Integer, List<T>>> mapCollection, final int numReferenceEdges) {
        final Map<Integer, List<T>> flattenedMap = new HashMap<>(SVUtils.hashMapCapacity(numReferenceEdges));
        for (final Map<Integer, List<T>> map : mapCollection) {
            for (final Map.Entry<Integer, List<T>> entry : map.entrySet()) {
                flattenedMap.putIfAbsent(entry.getKey(), new ArrayList<>());
                flattenedMap.get(entry.getKey()).addAll(entry.getValue());
            }
        }
        return flattenedMap;
    }

    private static void setGenotypeEvidenceProbabilities(final Collection<SVGraphGenotype> paths, final SVGraph graph) {
        final Collection<IndexedSVGraphEdge> breakpointEdges = graph.getEdges().stream().filter(e -> !e.isReference()).collect(Collectors.toList());
        for (final SVGraphGenotype genotype : paths) {
            final Map<Integer, Long> nonZeroEdgeVisits = genotype.getHaplotypes().stream()
                    .flatMap(h -> h.getEdges().stream())
                    .filter(e -> !e.isReference())
                    .collect(Collectors.groupingBy(edge -> edge.getIndex(), Collectors.counting()));
            final Set<Integer> zeroEdgeVisits = breakpointEdges.stream()
                    .map(IndexedSVGraphEdge::getIndex)
                    .filter(e -> !nonZeroEdgeVisits.containsKey(e))
                    .collect(Collectors.toSet());
            final Map<Integer, Long> edgeVisits = Stream.concat(nonZeroEdgeVisits.entrySet().stream().map(entry -> new Tuple2<>(entry.getKey(), entry.getValue())),
                zeroEdgeVisits.stream().map(e -> new Tuple2<>(e, Long.valueOf(0)))).collect(Collectors.toMap(pair -> pair._1, pair -> pair._2));
            final double logP = edgeVisits.entrySet().stream().mapToDouble(entry -> graph.getEdges().get(entry.getKey()).getPrior().getLogPrior((entry.getValue().intValue()))).sum();
            genotype.setEvidenceProbability(Math.exp(logP));
        }
        final double totalP = paths.stream().mapToDouble(SVGraphGenotype::getEvidenceProbability).sum();
        for (final SVGraphGenotype genotype : paths) {
            genotype.setEvidenceProbability(genotype.getEvidenceProbability() / totalP);
        }
    }

    /**
     * Integrates event probabilities over all haplotypes on each reference edge interval. Returns tuples of events and their probabilities.
     */
    private static Collection<Tuple2<CalledSVGraphEvent, Double>> integrateEdgeEvents(final Collection<SVGraphGenotype> paths, final SVGraph graph) {
        final IndexedSVGraphPath referencePath = new IndexedSVGraphPath(graph.getReferenceEdges());
        final Collection<Map<Integer, List<SVGraphEvent>>> edgeEventMapsCollection = paths.stream().map(path -> getHaplotypeEvents(path, graph, referencePath)).collect(Collectors.toList());
        final Map<Integer, List<SVGraphEvent>> edgeEventMap = flattenMaps(edgeEventMapsCollection, referencePath.size());
        final Collection<Tuple2<CalledSVGraphEvent, Double>> events = new ArrayList<>();
        for (final Map.Entry<Integer, List<SVGraphEvent>> entry : edgeEventMap.entrySet()) {
            final Map<CalledSVGraphEvent.Type, List<SVGraphEvent>> typeMap = entry.getValue().stream().collect(Collectors.groupingBy(SVGraphEvent::getType));
            for (final CalledSVGraphEvent.Type type : typeMap.keySet()) {
                final Collection<SVGraphEvent> typedEvents = typeMap.get(type);
                final double totalProbability = sumProbabilities(typedEvents);
                final SVGraphEvent firstEvent = typedEvents.iterator().next();
                final int groupId = firstEvent.getGroupId();
                final int pathId = firstEvent.getPathId();
                final CalledSVGraphEvent newEvent = new CalledSVGraphEvent(type, firstEvent.getInterval(), groupId, pathId, true, totalProbability);
                events.add(new Tuple2<>(newEvent, totalProbability));
            }
        }
        return events;
    }

    /**
     * Helper method for merging events
     */
    private static CalledSVGraphEvent mergeSortedEvents(final List<CalledSVGraphEvent> eventsToMerge) {
        if (eventsToMerge == null || eventsToMerge.isEmpty()) {
            throw new IllegalArgumentException("Events list cannot be empty or null");
        }
        final CalledSVGraphEvent firstEvent = eventsToMerge.get(0);
        final CalledSVGraphEvent lastEvent = eventsToMerge.get(eventsToMerge.size() - 1);
        final SVInterval firstInterval = firstEvent.getInterval();
        final SVInterval interval = new SVInterval(firstInterval.getContig(), firstInterval.getStart(), lastEvent.getInterval().getEnd());
        final CalledSVGraphEvent mergedEvent = new CalledSVGraphEvent(firstEvent.getType(), interval, firstEvent.getGroupId(), firstEvent.getPathId(), true, -1); //TODO set probability
        return mergedEvent;
    }

    /**
     * Finds and merges contiguous events
     */
    private static Collection<CalledSVGraphEvent> mergeAdjacentEvents(final Collection<CalledSVGraphEvent> events) {
        final Collection<CalledSVGraphEvent> mergedEvents = new ArrayList<>(events.size());
        final List<CalledSVGraphEvent> eventsToMerge = new ArrayList<>(events.size());
        final Set<CalledSVGraphEvent.Type> types = events.stream().map(CalledSVGraphEvent::getType).collect(Collectors.toSet());
        for (final CalledSVGraphEvent.Type type : types) {
            final List<CalledSVGraphEvent> sortedEventList = events.stream()
                    .filter(event -> event.getType().equals(type)).sorted((a, b) -> SVIntervalUtils.compareIntervals(a.getInterval(), b.getInterval()))
                    .collect(Collectors.toList());
            for (final CalledSVGraphEvent event : sortedEventList) {
                if (eventsToMerge.isEmpty()) {
                    eventsToMerge.add(event);
                } else {
                    final CalledSVGraphEvent previousEvent = eventsToMerge.get(eventsToMerge.size() - 1);
                    final SVInterval previousEventInterval = previousEvent.getInterval();
                    final SVInterval eventInterval = event.getInterval();
                    if (previousEventInterval.getContig() == eventInterval.getContig() &&
                            previousEventInterval.getEnd() == eventInterval.getStart()) {
                        eventsToMerge.add(event);
                    } else if (!eventInterval.equals(previousEventInterval)) {
                        mergedEvents.add(mergeSortedEvents(eventsToMerge));
                        eventsToMerge.clear();
                        eventsToMerge.add(event);
                    }
                }
            }
            if (!eventsToMerge.isEmpty()) {
                mergedEvents.add(mergeSortedEvents(eventsToMerge));
                eventsToMerge.clear();
            }
        }
        return mergedEvents;
    }

    /**
     * Creates event that could not be successfully resolved
     */
    private CalledSVGraphEvent getUnresolvedEvent(final SVInterval interval, final int groupId) {
        return new CalledSVGraphEvent(CalledSVGraphEvent.Type.UR, interval, groupId, 0, false, 0);
    }

    /**
     * Calls events over graph partitions
     */
    public Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> callEvents() {

        final SVGraphPartitioner graphPartitioner = new SVGraphPartitioner(graph);
        final BiFunction<SVInterval, SVInterval, Boolean> partitionFunction = (a, b) -> graphPartitioningFunction(a, b, arguments.paritionReciprocalOverlap);
        final List<SVGraph> graphPartitions = graphPartitioner.getIndependentSubgraphs(partitionFunction);

        final Collection<CalledSVGraphEvent> events = new ArrayList<>();
        final Collection<CalledSVGraphGenotype> haplotypes = new ArrayList<>();
        int partitionId = 0;
        for (final SVGraph partition : graphPartitions) {
            System.out.println("partition = " + partitionId + " / " + graphPartitions.size());
            final Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> result = generateEvents(partition, partitionId, arguments.minEventProb,
                    arguments.maxPathLengthFactor, arguments.maxEdgeVisits, copyNumberPosteriorsTree,
                    arguments.maxBranches, arguments.ploidy, arguments.minEventSize, arguments.minHaplotypeProb);
            if (result == null) {
                for (final SVInterval interval : partition.getContigIntervals()) {
                    events.add(getUnresolvedEvent(interval, partitionId));
                }
            } else {
                haplotypes.addAll(result._1);
                events.addAll(result._2);
            }
            partitionId++;
        }
        return new Tuple2<>(haplotypes, events);
    }

    /**
     * Contains event information including probability
     */
    private static final class SVGraphEvent {
        private final int groupId;
        private final int pathId;
        private final double probability;
        private final double evidenceProbability;
        private final CalledSVGraphEvent.Type type;
        private final SVInterval interval;
        private final boolean resolved; //False if solution not found

        public SVGraphEvent(final CalledSVGraphEvent.Type type, final SVInterval interval,
                            final int groupId, final int pathId, final double probability,
                            final double evidenceProbability, final boolean resolved) {
            Utils.nonNull(type, "Type cannot be null");
            Utils.nonNull(interval, "Interval cannot be null");
            this.type = type;
            this.interval = interval;
            this.groupId = groupId;
            this.pathId = pathId;
            this.probability = probability;
            this.evidenceProbability = evidenceProbability;
            this.resolved = resolved;
        }

        public double getProbability() {
            return probability;
        }

        public double getEvidenceProbability() {
            return evidenceProbability;
        }

        public boolean isResolved() {
            return resolved;
        }

        public int getGroupId() {
            return groupId;
        }

        public int getPathId() {
            return pathId;
        }

        public CalledSVGraphEvent.Type getType() {
            return type;
        }

        public SVInterval getInterval() {
            return interval;
        }
    }


    /**
     * Container for copy number posterior likelihoods
     */
    private final static class EdgeCopyNumberPosterior {
        private final double[] logPosteriors;

        public EdgeCopyNumberPosterior(final double[] logPosteriors) {
            this.logPosteriors = logPosteriors;
        }

        public double getLogPosterior(final int i) {
            if (i >= logPosteriors.length || i < 0) return Double.NEGATIVE_INFINITY;
            return logPosteriors[i];
        }

        public int numCopyNumberStates() {
            return logPosteriors.length;
        }
    }

}
