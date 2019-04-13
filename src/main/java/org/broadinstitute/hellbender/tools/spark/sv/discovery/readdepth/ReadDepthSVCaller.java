package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import org.broadinstitute.hellbender.tools.spark.sv.DiscoverVariantsFromReadDepthArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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

    private static boolean graphRepartitioningFunction(final SVInterval a, final SVInterval b, final double minReciprocalOverlap) {
        return SVIntervalUtils.hasReciprocalOverlap(a, b, minReciprocalOverlap);
    }

    /**
     * Main function that enumerates graph paths, integrates event probabilities, and produces calls
     */
    public static Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> generateEvents(final SVGraph graph, final int groupId, final double minEventProb,
                                                                                                           final double maxPathLengthFactor, final int maxEdgeVisits,
                                                                                                           final SVIntervalTree<SVCopyNumberInterval> copyNumberPosteriorsTree,
                                                                                                           final int maxQueueSize, final int baselineCopyNumber, final int minSize,
                                                                                                           final double minHaplotypeProb, final int maxBreakpointsPerHaplotype) {

        if (graph.getContigIntervals().iterator().next().overlaps(new SVInterval(2, 163045923, 163052669))) {
            int x = 0;
        }
        if (baselineCopyNumber == 0) return new Tuple2<>(Collections.emptyList(), Collections.emptyList());
        final SVGraphGenotyper searcher = new SVGraphGenotyper(graph);
        System.out.println("\tEnumerating haplotypes");
        final Collection<IndexedSVGraphPath> paths = searcher.enumerate(maxPathLengthFactor, maxEdgeVisits, maxQueueSize, maxBreakpointsPerHaplotype);
        if (paths == null) return null;
        System.out.println("\tEnumerating genotypes from " + paths.size() + " haplotypes");
        final List<SVGraphGenotype> genotypes = enumerateGenotypes(paths, graph, copyNumberPosteriorsTree, groupId, baselineCopyNumber);
        if (genotypes == null) return null;
        System.out.println("\tGetting ref genotype of " + genotypes.size() + " genotypes");

        SVGraphGenotype refGenotype = null;
        for (final SVGraphGenotype genotype : genotypes) {
            boolean isHomRef = true;
            for (final IndexedSVGraphPath haplotype : genotype.getHaplotypes()) {
                for (final IndexedSVGraphEdge edge : haplotype.getEdges()) {
                    if (!edge.isReference()) {
                        isHomRef = false;
                        break;
                    }
                }
            }
            if (isHomRef) {
                refGenotype = genotype;
                break;
            }
        }

        System.out.println("\tGetting max genotypes of " + genotypes.size() + " genotypes");
        final SVGraphGenotype finalRefGenotype = refGenotype; // For use in lambda
        final OptionalDouble optionalMaxP = genotypes.stream().filter(g -> g != finalRefGenotype).mapToDouble(SVGraphGenotype::getDepthLikelihood).max();
        final double maxP = optionalMaxP.orElse(refGenotype.getDepthLikelihood());
        final int MAX_GENOTYPES = 50;
        Collection<SVGraphGenotype> maxGenotypes = new ArrayList<>();
        for (int i = 0; i < genotypes.size(); i++) {
            final SVGraphGenotype genotype = genotypes.get(i);
            final double genotypeP = genotype.getDepthLikelihood();
            if (genotypeP == maxP && genotype != refGenotype) {
                maxGenotypes.add(genotype);
                if (maxGenotypes.size() == MAX_GENOTYPES) break;
            }
        }
        if (maxGenotypes.isEmpty()) maxGenotypes.add(refGenotype);
        final double maxPL = -10 * maxP/Math.log(10);
        final double refPL = -10 * refGenotype.getDepthLikelihood()/Math.log(10);
        final double quality = Math.max(Math.min(maxPL - refPL, 99), -99);

        System.out.println("\tGetting events of " + maxGenotypes.size() + " maximal genotypes");
        final Collection<CalledSVGraphEvent> events = new ArrayList(maxGenotypes.size());
        final double MIN_DEPTH_SUPPORT = 0.5;
        for (final SVGraphGenotype genotype : maxGenotypes) {
            final List<Collection<SVGraphEvent>> eventsList = getHaplotypeEvents(genotype, graph, new IndexedSVGraphPath(graph.getReferenceEdges()));
            for (final Collection<SVGraphEvent> haplotypeEvents : eventsList) {
                final Collection<CalledSVGraphEvent> calledEvents = haplotypeEvents.stream()
                        .map(event -> new CalledSVGraphEvent(event.getType(), event.getInterval(), event.getGroupId(), event.getPathId(), true, quality))
                        .collect(Collectors.toList());
                final Collection<CalledSVGraphEvent> mergedEvents = mergeAdjacentEvents(calledEvents);
                final Collection<CalledSVGraphEvent> sizeFilteredEvents = filterEventsBySize(mergedEvents, minSize);
                for (final CalledSVGraphEvent event : sizeFilteredEvents) {
                    Iterator<SVIntervalTree.Entry<SVCopyNumberInterval>> overlapIter = copyNumberPosteriorsTree.overlappers(event.getInterval());
                    int overlapLen = 0;
                    while (overlapIter.hasNext()) {
                        overlapLen += overlapIter.next().getInterval().getLength();
                    }
                    final double overlapFraction = overlapLen / (double) event.getInterval().getLength();
                    if (overlapFraction >= MIN_DEPTH_SUPPORT) {
                        events.add(event);
                    }
                }
            }
        }

        if (graph.getContigIntervals().iterator().next().overlaps(new SVInterval(2, 163045923, 163052669))) {
            int x = 0;
        }

        final Collection<CalledSVGraphGenotype> haplotypes = convertToCalledHaplotypes(maxGenotypes, graph);
        System.out.println("\tReturning " + haplotypes.size() + " haplotypes and " + events.size() + " events");
        return new Tuple2<>(haplotypes, events);
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

    private static List<SVGraphGenotype> enumerateGenotypes(final Collection<IndexedSVGraphPath> paths,
                                                                  final SVGraph graph,
                                                                  final SVIntervalTree<SVCopyNumberInterval> copyNumberPosteriorsTree,
                                                                  final int groupId,
                                                                  final int baselineCopyNumber) {
        final int numEdges = graph.getEdges().size();
        final List<EdgeCopyNumberPosterior> copyNumberPosteriors = getEdgeCopyNumberPosteriors(copyNumberPosteriorsTree, graph);

        //If large number of copies, only allow a handful to be non-reference
        final int ignoredRefCopies = Math.max(baselineCopyNumber - 4, 0);
        final int nonRefCopies = Math.min(baselineCopyNumber, 4);
        if (4e9 < Math.pow(paths.size(), nonRefCopies)) {
            System.out.println("\tToo many genotype combinations: " + paths.size() + " ^ " + nonRefCopies);
            return null; //Genotypes wouldn't fit into an array
        }
        final List<IndexedSVGraphPath> pathsList = new ArrayList<>(paths);
        final List<SVGraphGenotype> genotypes = new ArrayList<>(paths.size() * paths.size());
        final List<int[]> edgeCopyNumberStates = pathsList.stream().map(path -> getEdgeCopyNumberStates(path, numEdges)).collect(Collectors.toList());

        //Enumerate index combinations (with repetition)
        final List<List<Integer>> combinations = new ArrayList<>();
        combinationsWithRepetition(0, paths.size()-1, nonRefCopies, new ArrayList<>(paths.size()), combinations);

        // Calculate posterior for each combination
        int genotypeId = 0;
        for (final List<Integer> combination : combinations) {
            final List<IndexedSVGraphPath> genotypePaths = new ArrayList<>(combination.size());
            final List<int[]> haplotypeStates = new ArrayList<>(combination.size());
            for (final Integer i : combination) {
                genotypePaths.add(pathsList.get(i));
                haplotypeStates.add(edgeCopyNumberStates.get(i));
            }
            final int[] genotypeStates = new int[numEdges];
            for (int i = 0; i < haplotypeStates.size(); i++) {
                for (int k = 0; k < numEdges; k++) {
                    genotypeStates[k] += haplotypeStates.get(i)[k];
                }
            }
            if (ignoredRefCopies > 0) {
                for (int k = 0; k < numEdges; k++) {
                    genotypeStates[k] += ignoredRefCopies;
                }
            }
            double logLikelihood = 0;
            for (int k = 0; k < copyNumberPosteriors.size(); k++) {
                logLikelihood += copyNumberPosteriors.get(k).getLogPosterior(genotypeStates[k]);
            }

            final SVGraphGenotype genotype = new SVGraphGenotype(groupId, genotypeId, genotypePaths);
            genotype.setDepthLikelihood(logLikelihood);
            genotypes.add(genotype);
            genotypeId++;
        }
        return genotypes;
    }

    static void combinationsWithRepetition(int start, int end, int r, List<Integer> currentCombination, List<List<Integer>> combinations) {
        if (r == 0) {
            combinations.add(new ArrayList<>(currentCombination));
        } else {
            for (int i = start; i <= end; i++) {
                currentCombination.add(i);
                combinationsWithRepetition(i, end, r - 1, currentCombination, combinations);
                currentCombination.remove(currentCombination.size() - 1);
            }
        }
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
            final double p = denom == 0 ? 0 : genotype.getEvidenceProbability() * genotype.getDepthProbability() / denom;
            genotype.setProbability(p);
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
        boolean deletion = false;
        for (int i = 0; i < referenceEdgeCountsList.size(); i++) {
            if (referenceEdgeCountsList.get(i)[edgeIndex] < 1) {
                deletion = true;
                break;
            }
        }
        boolean duplication = false;
        for (int i = 0; i < referenceEdgeCountsList.size(); i++) {
            if (referenceEdgeCountsList.get(i)[edgeIndex] > 1) {
                duplication = true;
                break;
            }
        }
        boolean inversion = false;
        for (int i = 0; i < referenceEdgeInversionsList.size(); i++) {
            if (referenceEdgeInversionsList.get(i)[edgeIndex]) {
                inversion = true;
                break;
            }
        }
        final int numEvents = (deletion ? 1 : 0) + (duplication ? 1 : 0) + (inversion ? 1 : 0) - (deletion && duplication ? 1 : 0);
        if (numEvents == 0) return Collections.emptyList();
        final List<SVGraphEvent> events = new ArrayList<>(numEvents);
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
        return events;
    }

    /**
     * Aggregates potential event calls for the given haplotypes, producing a map from reference edge index to list of events
     */
    private static List<Collection<SVGraphEvent>> getHaplotypeEvents(final SVGraphGenotype haplotypes, final SVGraph graph, final IndexedSVGraphPath referencePath) {
        //System.out.println("\t\t\tGetting events for haplotypes " + haplotypes.getGenotypeId());
        final List<IndexedSVGraphEdge> edges = graph.getEdges();
        //System.out.println("\t\t\t\tgetReferenceEdgeCountsAndInversions");
        final Tuple2<List<int[]>, List<boolean[]>> referenceEdgeResults = getReferenceEdgeCountsAndInversions(haplotypes, edges);
        final List<int[]> referenceEdgeCountsList = referenceEdgeResults._1;
        final List<boolean[]> referenceEdgeInversionsList = referenceEdgeResults._2;

        final List<SVGraphNode> nodes = graph.getNodes();
        final List<Collection<SVGraphEvent>> events = new ArrayList<>(haplotypes.getHaplotypes().size());
        final int groupId = haplotypes.getGroupId();
        final int pathId = haplotypes.getGenotypeId();
        final double probability = haplotypes.getDepthProbability();
        final double evidenceProbability = haplotypes.getEvidenceProbability();
        //System.out.println("\t\t\t\tgetEdgeEvents loop over " + referencePath.getEdges().size() + " edges");
        for (int i = 0; i < haplotypes.getHaplotypes().size(); i++) {
            final int[] edgeCounts = referenceEdgeCountsList.get(i);
            final boolean[] inversionsList = referenceEdgeInversionsList.get(i);
            final List<SVGraphEvent> haplotypeEvents = new ArrayList<>();
            for (int j = 0; j < edgeCounts.length; j++) {
                if (edges.get(j).isReference()) {
                    if (edgeCounts[j] == 0) {
                        haplotypeEvents.add(new SVGraphEvent(CalledSVGraphEvent.Type.DEL, edges.get(j).getInterval(), groupId, pathId, probability, evidenceProbability, true));
                    } else if (edgeCounts[j] > 1) {
                        haplotypeEvents.add(new SVGraphEvent(CalledSVGraphEvent.Type.DUP, edges.get(j).getInterval(), groupId, pathId, probability, evidenceProbability, true));
                    }
                    if (inversionsList[j]) {
                        haplotypeEvents.add(new SVGraphEvent(CalledSVGraphEvent.Type.INV, edges.get(j).getInterval(), groupId, pathId, probability, evidenceProbability, true));
                    }
                }
            }
            events.add(haplotypeEvents);
        }
        return events;
    }

    private static double sumProbabilities(final Collection<SVGraphEvent> events) {
        return events.stream().mapToDouble(SVGraphEvent::getProbability).sum();
    }

    private static <T> Collection<Collection<T>> flattenLists(final Collection<List<Collection<T>>> mapCollection) {
        if (mapCollection.isEmpty()) return Collections.emptyList();
        final int maxSize = mapCollection.stream().mapToInt(list -> list.size()).max().getAsInt();
        final List<Collection<T>> flattenedList = new ArrayList<>(maxSize);
        for (int i = 0; i < maxSize; i++) {
            flattenedList.add(new ArrayList<>());
        }
        for (final List<Collection<T>> entry : mapCollection) {
            for (int i = 0; i < entry.size(); i++) {
                if (entry.get(i) != null) {
                    flattenedList.get(i).addAll(entry.get(i));
                }
            }
        }
        return flattenedList;
    }

    private static void setGenotypeEvidenceProbabilities(final Collection<SVGraphGenotype> paths, final SVGraph graph) {
        /*
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
        */
        for (final SVGraphGenotype genotype : paths) {
            genotype.setEvidenceProbability(1);
        }
    }

    /**
     * Integrates event probabilities over all haplotypes on each reference edge interval. Returns tuples of events and their probabilities.
     */
    private static Collection<Tuple2<CalledSVGraphEvent, Double>> integrateEdgeEvents(final Collection<SVGraphGenotype> paths, final SVGraph graph) {
        final IndexedSVGraphPath referencePath = new IndexedSVGraphPath(graph.getReferenceEdges());
        //System.out.println("\t\tCollecting haplotype events...");
        final Collection<List<Collection<SVGraphEvent>>> edgeEventMapsCollection = paths.stream().map(path -> getHaplotypeEvents(path, graph, referencePath)).collect(Collectors.toList());
        //System.out.println("\t\tFlattening maps...");
        final Collection<Collection<SVGraphEvent>> edgeEvents = flattenLists(edgeEventMapsCollection);
        final Collection<Tuple2<CalledSVGraphEvent, Double>> events = new ArrayList<>();
        int i = 0;
        for (final Collection<SVGraphEvent> value : edgeEvents) {
            //System.out.println("\t\tIntegrating entry " + ++i + " / " + edgeEvents.size());
            final Map<CalledSVGraphEvent.Type, List<SVGraphEvent>> typeMap = value.stream().collect(Collectors.groupingBy(SVGraphEvent::getType));
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
        final CalledSVGraphEvent mergedEvent = new CalledSVGraphEvent(firstEvent.getType(), interval, firstEvent.getGroupId(), firstEvent.getPathId(), true, firstEvent.getProbability());
        return mergedEvent;
    }

    /**
     * Finds and merges contiguous events of equal probability and type
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
                            previousEventInterval.getEnd() == eventInterval.getStart() &&
                            previousEvent.getProbability() == event.getProbability()) {
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

    private List<List<Integer>> getPartitionDependenceGraph(final List<SVGraph> graphPartitions) {
        //Create interval tree of the partitions
        final SVIntervalTree<Integer> graphTree = new SVIntervalTree<>();
        for (int i = 0; i < graphPartitions.size(); i++) {
            for (final SVInterval interval : graphPartitions.get(i).getContigIntervals()) {
                graphTree.put(interval, i);
            }
        }

        //Create dependence graph matrix, in which i -> j if i is the smallest partition that completely overlaps (and is larger than) j
        final List<List<Integer>> dependenceGraph = new ArrayList<>(graphPartitions.size());
        final List<List<Integer>> parentGraph = new ArrayList<>(graphPartitions.size());
        for (int i = 0; i < graphPartitions.size(); i++) {
            dependenceGraph.add(new ArrayList<>());
            parentGraph.add(new ArrayList<>());
        }
        for (int i = 0; i < graphPartitions.size(); i++) {
            final Collection<SVInterval> iIntervals = graphPartitions.get(i).getContigIntervals();
            for (final SVInterval iInterval : iIntervals) {
                // Find smallest parent on this contig
                final Iterator<SVIntervalTree.Entry<Integer>> iter = graphTree.overlappers(iInterval);
                int minSize = Integer.MAX_VALUE;
                int minParent = -1;
                while (iter.hasNext()) {
                    final SVIntervalTree.Entry<Integer> entry = iter.next();
                    final int j = entry.getValue();
                    if (j == i) continue;    //Don't introduce self-edges
                    final SVInterval jInterval = entry.getInterval();
                    if (jInterval.getLength() < minSize) {
                        minParent = j;
                    }
                }
                if (minParent != -1) {
                    dependenceGraph.get(minParent).add(i);
                    parentGraph.get(i).add(minParent);
                }
            }
        }

        return dependenceGraph.stream().map(ArrayList::new).collect(Collectors.toList());
    }

    /**
     * Calls events over graph partitions
     */
    public Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> callEvents() {
        final BiFunction<SVInterval, SVInterval, Boolean> partitionFunction = (a, b) -> graphPartitioningFunction(a, b, arguments.paritionReciprocalOverlap);
        final SVGraphPartitioner graphPartitioner = new SVGraphPartitioner(graph);

        // Partition the graph
        final List<SVGraph> graphPartitions = graphPartitioner.getIndependentSubgraphs(partitionFunction);

        // Determine partition hierarchy (partitions that contain others)
        final List<List<Integer>> partitionDependenceGraph = getPartitionDependenceGraph(graphPartitions);
        final Set<Integer> childPartitions = partitionDependenceGraph.stream().flatMap(List::stream).collect(Collectors.toSet());
        final List<Integer> rootPartitions = IntStream.range(0, graphPartitions.size()).filter(i -> childPartitions.contains(i)).boxed().collect(Collectors.toList());

        // Traverse the graph starting at each root
        final int numPartitions = graphPartitions.size();
        final boolean[] processedPartitions = new boolean[numPartitions];
        final Collection<CalledSVGraphGenotype> haplotypes = new ArrayList<>();
        final Collection<CalledSVGraphEvent> events = new ArrayList<>();
        for (final Integer root : rootPartitions) {
            final Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> result = callEventsHelper(root, graphPartitions, processedPartitions, Collections.emptyList(), 2, partitionDependenceGraph);
            haplotypes.addAll(result._1);
            events.addAll(result._2);
        }
        return new Tuple2<>(haplotypes, events);
    }

    //TODO : does not process inter-chromosomal partitions properly
    private Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> callEventsHelper(final int partitionIndex,
                                                                                                       final List<SVGraph> partitions,
                                                                                                       final boolean[] processedPartitions,
                                                                                                       final Collection<CalledSVGraphGenotype> parentHaplotypes,
                                                                                                       final int defaultCopyNumber,
                                                                                                       final List<List<Integer>> partitionDependenceGraph) {
        if (processedPartitions[partitionIndex]) {
            return new Tuple2<>(Collections.emptyList(), Collections.emptyList());
        }
        processedPartitions[partitionIndex] = true;

        final Collection<CalledSVGraphGenotype> haplotypes = new ArrayList<>();
        final Collection<CalledSVGraphEvent> events = new ArrayList<>();
        final int baselineCopyNumber = defaultCopyNumber;
        final List<CalledSVGraphGenotype> partitionHaplotypes = new ArrayList<>();
        final SVGraph partition = partitions.get(partitionIndex);

        final Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> result = generateEvents(partition, partitionIndex,
                arguments.minEventProb, arguments.maxPathLengthFactor, arguments.maxEdgeVisits, copyNumberPosteriorsTree,
                arguments.maxBranches, baselineCopyNumber, arguments.minEventSize, arguments.minHaplotypeProb, Integer.MAX_VALUE);
        if (result == null) {
            final BiFunction<SVInterval, SVInterval, Boolean> repartitionFunction = (a, b) -> graphRepartitioningFunction(a, b, 0.1);
            final SVGraphPartitioner repartitioner = new SVGraphPartitioner(partition);
            final List<SVGraph> repartitions = repartitioner.getIndependentSubgraphs(repartitionFunction);
            for (final SVGraph repartition : repartitions) {
                final Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> retryResult = generateEvents(repartition, partitionIndex,
                        arguments.minEventProb, arguments.maxPathLengthFactor, arguments.maxEdgeVisits, copyNumberPosteriorsTree,
                        arguments.maxBranches, baselineCopyNumber, arguments.minEventSize, arguments.minHaplotypeProb, Integer.MAX_VALUE);
                if (retryResult == null) {
                    for (final IndexedSVGraphEdge edge : repartition.getEdges()) {
                        if (!edge.isReference()) {
                            events.add(getUnresolvedEvent(edge.getInterval(), partitionIndex));
                        }
                    }
                } else {
                    partitionHaplotypes.addAll(retryResult._1);
                    haplotypes.addAll(retryResult._1);
                    events.addAll(retryResult._2);
                }
            }
        } else {
            partitionHaplotypes.addAll(result._1);
            haplotypes.addAll(result._1);
            events.addAll(result._2);
        }

        for (final Integer childIndex : partitionDependenceGraph.get(partitionIndex)) {
            //TODO find correct baseline copy state
            final Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> childResult = callEventsHelper(childIndex, partitions, processedPartitions, partitionHaplotypes, baselineCopyNumber, partitionDependenceGraph);
            haplotypes.addAll(childResult._1);
            events.addAll(childResult._2);
        }
        return new Tuple2<>(haplotypes, events);

            /*System.out.println("\tSecond attempt failed on " + partitionInterval + " (N = " + partition.getNodes().size() + ", E = " + partition.getEdges().size() + ", E_ref = " + partition.getReferenceEdges().size() + ")");
            for (final IndexedSVGraphEdge edge : partition.getEdges()) {
                if (!edge.isReference()) {
                    events.add(getUnresolvedEvent(edge.getInterval(), partitionIndex));
                }
            }*/


        /*
        if (!parentHaplotypes.isEmpty()) {
            int maxProbHaplotypeIndex = -1;
            double maxHaplotypeProb = -1;
            for (int i = 0; i < parentHaplotypes.size(); i++) {
                final double prob = parentHaplotypes.get(i).getProbability();
                if (prob > maxHaplotypeProb) {
                    maxHaplotypeProb = prob;
                    maxProbHaplotypeIndex = i;
                }
            }
            final CalledSVGraphGenotype maxProbHaplotype = parentHaplotypes.get(maxProbHaplotypeIndex);
            baselineCopyNumber = (int) maxProbHaplotype.getHaplotypes().stream()
                    .flatMap(haplotype -> haplotype.getEdges().stream().filter(edge -> edge.isReference() && edge.getInterval().overlaps(partitionInterval)))
                    .count();
        } else {
            baselineCopyNumber = defaultCopyNumber;
        }
        */
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
