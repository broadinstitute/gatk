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
                                                                                                           final double minHaplotypeProb, final int maxBreakpointsPerHaplotype,
                                                                                                           final Double parentQuality) {

        if (SVIntervalUtils.hasReciprocalOverlap(graph.getContigIntervals().iterator().next(), new SVInterval(6, 113776106, 113782152), 0.5)) {
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
                if (!haplotype.isReference()) {
                    isHomRef = false;
                    break;
                }
            }
            if (isHomRef) {
                refGenotype = genotype;
                break;
            }
        }

        //System.out.println("\tGetting max genotypes of " + genotypes.size() + " genotypes");
        final int MAX_GENOTYPES = 100000;
        final double refP = refGenotype.getDepthLikelihood();

        /*
        final SVGraphGenotype finalRefGenotype = refGenotype; // For use in lambda
        final OptionalDouble optionalMaxP = genotypes.stream().filter(g -> g != finalRefGenotype).mapToDouble(SVGraphGenotype::getDepthLikelihood).max();
        final double maxP = optionalMaxP.orElse(refGenotype.getDepthLikelihood());
        List<SVGraphGenotype> maxGenotypes = new ArrayList<>();
        for (int i = 0; i < genotypes.size(); i++) {
            final SVGraphGenotype genotype = genotypes.get(i);
            final double genotypeP = genotype.getDepthLikelihood();
            if ((genotypeP == maxP || genotypeP >= refP) && genotype != refGenotype) {
                maxGenotypes.add(genotype);
                if (maxGenotypes.size() == MAX_GENOTYPES) break;
            }
        }
        if (maxGenotypes.isEmpty()) maxGenotypes.add(refGenotype);
        */

        if (genotypes.size() > MAX_GENOTYPES) return null; //TODO give up if too many genotypes
        final int numGenotypes = Math.min(genotypes.size(), MAX_GENOTYPES);
        System.out.println("\tGetting events of " + numGenotypes + " maximal genotypes");
        final List<IndexedSVGraphEdge> referenceEdges = graph.getReferenceEdges();
        final Collection<CalledSVGraphEvent> events = new ArrayList(numGenotypes);
        final double MIN_DEPTH_SUPPORT = 0.5;
        final double refPL = -10 * refP/Math.log(10);
        Collections.sort(genotypes, Comparator.comparingDouble(g -> -g.getDepthLikelihood()));
        for (int i = 0; i < numGenotypes; i++) {
            final SVGraphGenotype genotype = genotypes.get(i);
            final double genotypePL = -10 * genotype.getDepthLikelihood()/Math.log(10);
            final double genotypeQuality = Math.max(Math.min(genotypePL - refPL, 99), -99);
            final double quality = parentQuality != null ? Math.max(genotypeQuality, parentQuality) : genotypeQuality;
            genotype.setProbability(quality);
            if (quality >= 99) continue;
            final Collection<SVGraphEvent> eventsCollection = getHaplotypeEvents(genotype, graph, new IndexedSVGraphPath(referenceEdges));
            final Collection<CalledSVGraphEvent> calledEvents = new ArrayList<>(eventsCollection.size());
            for (final SVGraphEvent event : eventsCollection) {
                calledEvents.add(new CalledSVGraphEvent(event.getType(), event.getInterval(), event.getGroupId(), event.getPathId(), true, quality, event.isHomozygous()));
            }
            final Collection<CalledSVGraphEvent> mergedEvents = mergeAdjacentEvents(calledEvents);
            final Collection<CalledSVGraphEvent> sizeFilteredEvents = filterEventsBySize(mergedEvents, minSize);
            events.addAll(sizeFilteredEvents);
            /*
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
            }*/
        }

        // Deduplicate events, taking the most likely one
        // TODO treats hom/het the same
        System.out.println("\tDeduplicating " + events.size() + " events.");
        final Map<CalledSVGraphEvent.Type,Map<SVInterval,List<CalledSVGraphEvent>>> eventsMap = new HashMap<>(SVUtils.hashMapCapacity(CalledSVGraphEvent.Type.values().length));
        for (final CalledSVGraphEvent.Type type : CalledSVGraphEvent.Type.values()) {
            eventsMap.put(type, new HashMap<>());
        }
        for (final CalledSVGraphEvent event : events) {
            eventsMap.get(event.getType()).putIfAbsent(event.getInterval(), new ArrayList<>());
            eventsMap.get(event.getType()).get(event.getInterval()).add(event);
        }
        final Collection<CalledSVGraphEvent> collapsedEvents = new ArrayList<>(eventsMap.size());
        for (final CalledSVGraphEvent.Type type : eventsMap.keySet()){
            for (final List<CalledSVGraphEvent> list : eventsMap.get(type).values()) {
                double minQuality = 99;
                for (final CalledSVGraphEvent event : list) {
                    if (event.getRefQuality() < minQuality) minQuality = event.getRefQuality();
                }
                for (final CalledSVGraphEvent event : list) {
                    if (event.getRefQuality() == minQuality) {
                        collapsedEvents.add(event);
                        break;
                    }
                }
            }
        }

        //Compute GQs
        /*
        for (final CalledSVGraphEvent event : collapsedEvents) {
            double total = 0;
            for (final CalledSVGraphEvent otherEvent : collapsedEvents) {
                if (otherEvent != event && otherEvent.getPathId() != event.getPathId()) {
                    final double overlap = event.getInterval().overlapLen(otherEvent.getInterval()) / (double) event.getInterval().getLength();
                    total += overlap * Math.max(0, event.getRefQuality() - otherEvent.getRefQuality());
                }
            }
            event.setGenotypeQuality(Math.min(total, 99));
        }*/

        //TODO this could be better - events should be organized by genotype
        for (final CalledSVGraphEvent event : collapsedEvents) {
            double total = 0;
            for (final CalledSVGraphEvent otherEvent : collapsedEvents) {
                if (otherEvent.getPathId() != event.getPathId()) {
                    final double overlap = event.getInterval().overlapLen(otherEvent.getInterval()) / (double) event.getInterval().getLength();
                    total += overlap * Math.exp(-otherEvent.getRefQuality() / 10.0);
                }
            }
            final double eventTerm = Math.exp(-event.getRefQuality() / 10.0); //Add at the end so we don't penalize for other events in the genotype
            final double p = eventTerm / (total + eventTerm);
            final double gq = -10.0 * Math.log10(Math.max(p, Double.MIN_VALUE));
            event.setGenotypeQuality(Math.min(gq, 99));
        }

        final Collection<CalledSVGraphGenotype> haplotypes = convertToCalledHaplotypes(genotypes, graph);
        System.out.println("\tReturning " + haplotypes.size() + " haplotypes and " + collapsedEvents.size() + " events");
        return new Tuple2<>(haplotypes, collapsedEvents);
    }

    private static Collection<CalledSVGraphGenotype> convertToCalledHaplotypes(final Collection<SVGraphGenotype> haplotypes, final SVGraph graph) {
        return haplotypes.stream().map(h -> new CalledSVGraphGenotype(h, graph)).collect(Collectors.toList());
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
        final int ignoredRefCopies = Math.max(baselineCopyNumber - 2, 0);
        final int nonRefCopies = Math.min(baselineCopyNumber, 2);
        if (4e9 < Math.pow(paths.size(), nonRefCopies)) {
            System.out.println("\tToo many genotype combinations: " + paths.size() + " ^ " + nonRefCopies);
            return null; //Genotypes wouldn't fit into an array
        }

        final List<IndexedSVGraphPath> pathsList = new ArrayList<>(paths);
        int refPath = -1;
        for (int i = 0; i < pathsList.size(); i++) {
            if (pathsList.get(i).isReference()) {
                refPath = i;
                break;
            }
        }
        //TODO check if refpath is -1

        final List<int[]> edgeCopyNumberStates = pathsList.stream().map(path -> getEdgeCopyNumberStates(path, numEdges)).collect(Collectors.toList());
        final int PATH_LIMIT = 1000;
        final int numGenotypes = paths.size() < PATH_LIMIT ? paths.size() * paths.size() : paths.size() * nonRefCopies;
        final List<SVGraphGenotype> genotypes = new ArrayList<>(numGenotypes);

        //Enumerate index combinations (with repetition)
        final List<List<Integer>> combinations = new ArrayList<>();
        combinationsWithRepetition(0, paths.size() - 1, nonRefCopies, new ArrayList<>(paths.size()), combinations);

        // Calculate posterior for each combination
        int genotypeId = 0;
        for (final List<Integer> combination : combinations) {
            final Set<Integer> nonRefAlleles = new HashSet<>(SVUtils.hashMapCapacity(combination.size()));
            for (final Integer index : combination) {
                if (index != refPath) {
                    nonRefAlleles.add(index);
                }
            }
            //TODO only check single-allele combinations if there are too many haplotypes
            if (paths.size() < PATH_LIMIT || nonRefAlleles.size() <= 1) {
                final List<IndexedSVGraphPath> genotypePaths = new ArrayList<>(combination.size());
                final List<int[]> haplotypeStates = new ArrayList<>(combination.size());
                for (final Integer i : combination) {
                    genotypePaths.add(pathsList.get(i));
                    haplotypeStates.add(edgeCopyNumberStates.get(i));
                }
                final SVGraphGenotype genotype = getGenotypeCombination(haplotypeStates, copyNumberPosteriors, genotypePaths, numEdges, ignoredRefCopies, groupId, genotypeId);
                genotypes.add(genotype);
                genotypeId++;
            }
        }
        return genotypes;
    }

    private static SVGraphGenotype getGenotypeCombination(final List<int[]> haplotypeStates,
                                                          final List<EdgeCopyNumberPosterior> copyNumberPosteriors,
                                                          final List<IndexedSVGraphPath> genotypePaths,
                                                          final int numEdges, final int ignoredRefCopies,
                                                          final int groupId, final int genotypeId) {
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
        return genotype;
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
     * Aggregates potential event calls for the given haplotypes, producing a map from reference edge index to list of events
     */
    private static Collection<SVGraphEvent> getHaplotypeEvents(final SVGraphGenotype haplotypes, final SVGraph graph, final IndexedSVGraphPath referencePath) {
        //System.out.println("\t\t\tGetting events for haplotypes " + haplotypes.getGenotypeId());
        final List<IndexedSVGraphEdge> edges = graph.getEdges();
        //System.out.println("\t\t\t\tgetReferenceEdgeCountsAndInversions");
        final Tuple2<List<int[]>, List<boolean[]>> referenceEdgeResults = getReferenceEdgeCountsAndInversions(haplotypes, edges);
        final List<int[]> referenceEdgeCountsList = referenceEdgeResults._1;
        final List<boolean[]> referenceEdgeInversionsList = referenceEdgeResults._2;

        final Collection<SVGraphEvent> rawEvents = new ArrayList<>(haplotypes.getHaplotypes().size());
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
                        haplotypeEvents.add(new SVGraphEvent(CalledSVGraphEvent.Type.DEL, edges.get(j).getInterval(), groupId, pathId, probability, evidenceProbability, true, false));
                    } else if (edgeCounts[j] > 1) {
                        haplotypeEvents.add(new SVGraphEvent(CalledSVGraphEvent.Type.DUP, edges.get(j).getInterval(), groupId, pathId, probability, evidenceProbability, true, false));
                    }
                    if (inversionsList[j]) {
                        haplotypeEvents.add(new SVGraphEvent(CalledSVGraphEvent.Type.INV, edges.get(j).getInterval(), groupId, pathId, probability, evidenceProbability, true, false));
                    }
                }
            }
            rawEvents.addAll(haplotypeEvents);
        }

        // Find homozygous events
        final Map<CalledSVGraphEvent.Type,Map<SVInterval,List<SVGraphEvent>>> eventsMap = new HashMap<>(SVUtils.hashMapCapacity(CalledSVGraphEvent.Type.values().length));
        for (final CalledSVGraphEvent.Type type : CalledSVGraphEvent.Type.values()) {
            eventsMap.put(type, new HashMap<>());
        }
        for (final SVGraphEvent event : rawEvents) {
            eventsMap.get(event.getType()).putIfAbsent(event.getInterval(), new ArrayList<>());
            eventsMap.get(event.getType()).get(event.getInterval()).add(event);
        }
        final List<SVGraphEvent> genotypedEvents = new ArrayList<>(rawEvents.size());
        for (final CalledSVGraphEvent.Type type : eventsMap.keySet()) {
            for (final List<SVGraphEvent> eventList : eventsMap.get(type).values()) {
                final SVGraphEvent event = eventList.get(0);
                if (eventList.size() == 1) {
                    genotypedEvents.add(event);
                } else {
                    genotypedEvents.add(new SVGraphEvent(type, event.getInterval(), event.getGroupId(), event.getPathId(), event.getProbability(), event.getEvidenceProbability(), event.isResolved(), true));
                }
            }
        }
        return genotypedEvents;
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
        final CalledSVGraphEvent mergedEvent = new CalledSVGraphEvent(firstEvent.getType(), interval, firstEvent.getGroupId(), firstEvent.getPathId(), true, firstEvent.getRefQuality(), firstEvent.isHomozygous());
        return mergedEvent;
    }

    /**
     * Finds and merges contiguous events of equal probability, type, and zygosity
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
                            previousEvent.getRefQuality() == event.getRefQuality() &&
                            previousEvent.isHomozygous() == event.isHomozygous()) {
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
        return new CalledSVGraphEvent(CalledSVGraphEvent.Type.UR, interval, groupId, 0, false, 0, false);
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
                    final SVInterval jInterval = entry.getInterval();
                    if (iInterval.getLength() >= jInterval.getLength()) continue; // Make sure i is strictly smaller
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

        for (int i = 0; i < graphPartitions.size(); i++) {
            for (final IndexedSVGraphEdge edge : graphPartitions.get(i).getEdges()) {
                if (SVIntervalUtils.hasReciprocalOverlap(edge.getInterval(), new SVInterval(6, 113776106, 113782152), 0.1)) {
                    int x = 0;
                }
            }
        }


        // Determine partition hierarchy (partitions that contain others)
        final List<List<Integer>> partitionDependenceGraph = getPartitionDependenceGraph(graphPartitions);
        final Set<Integer> childPartitions = partitionDependenceGraph.stream().flatMap(List::stream).collect(Collectors.toSet());
        final List<Integer> rootPartitions = IntStream.range(0, graphPartitions.size()).filter(i -> !childPartitions.contains(i)).boxed().collect(Collectors.toList());

        // Traverse the graph starting at each root
        final int numPartitions = graphPartitions.size();
        final boolean[] processedPartitions = new boolean[numPartitions];
        final Collection<CalledSVGraphGenotype> haplotypes = new ArrayList<>();
        final Collection<CalledSVGraphEvent> events = new ArrayList<>();
        for (final Integer root : rootPartitions) {
            final Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> result = callEventsHelper(root, graphPartitions, processedPartitions, arguments.ploidy, null, partitionDependenceGraph);
            haplotypes.addAll(result._1);
            events.addAll(result._2);
        }
        return new Tuple2<>(haplotypes, events);
    }

    //TODO : does not process inter-chromosomal partitions properly
    private Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> callEventsHelper(final int partitionIndex,
                                                                                                       final List<SVGraph> partitions,
                                                                                                       final boolean[] processedPartitions,
                                                                                                       final int baselineCopyNumber,
                                                                                                       final Double parentQuality,
                                                                                                       final List<List<Integer>> partitionDependenceGraph) {

        final SVGraph partition = partitions.get(partitionIndex);

        for (final IndexedSVGraphEdge edge : partition.getEdges()) {
            if (SVIntervalUtils.hasReciprocalOverlap(edge.getInterval(), new SVInterval(6, 113776106, 113782152), 0.1)) {
                int x = 0;
            }
        }

        if (processedPartitions[partitionIndex]) {
            return new Tuple2<>(Collections.emptyList(), Collections.emptyList());
        }
        processedPartitions[partitionIndex] = true;

        //System.out.println("Partition " + partitionIndex + " / " + partitions.size());

        final Collection<CalledSVGraphGenotype> haplotypes = new ArrayList<>();
        final Collection<CalledSVGraphEvent> events = new ArrayList<>();
        final List<CalledSVGraphGenotype> partitionHaplotypes = new ArrayList<>();

        final Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> result = generateEvents(partition, partitionIndex,
                arguments.minEventProb, arguments.maxPathLengthFactor, arguments.maxEdgeVisits, copyNumberPosteriorsTree,
                arguments.maxBranches, baselineCopyNumber, arguments.minEventSize, arguments.minHaplotypeProb, Integer.MAX_VALUE, parentQuality);
        if (result == null) {
            /*System.out.println("\tUnresolved");
            for (final IndexedSVGraphEdge edge : partition.getEdges()) {
                if (!edge.isReference()) {
                    events.add(getUnresolvedEvent(edge.getInterval(), partitionIndex));
                }
            }*/

            for (int maxBreakpoints = 2; maxBreakpoints >= 1; maxBreakpoints -= 1) {
                System.out.println("Retrying with maxBreakpoints = " + maxBreakpoints);
                final Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> retryResult = generateEvents(partition, partitionIndex,
                        arguments.minEventProb, arguments.maxPathLengthFactor, arguments.maxEdgeVisits, copyNumberPosteriorsTree,
                        arguments.maxBranches, baselineCopyNumber, arguments.minEventSize, arguments.minHaplotypeProb, maxBreakpoints, parentQuality);
                if (retryResult != null) {
                    System.out.println("Retry success at maxBreakpoints = " + maxBreakpoints);
                    partitionHaplotypes.addAll(retryResult._1);
                    //TODO this list gets too big
                    //haplotypes.addAll(retryResult._1);
                    events.addAll(retryResult._2);
                    break;
                } else if (maxBreakpoints == 1) {
                    System.out.println("\tUnresolved");
                    for (final IndexedSVGraphEdge edge : partition.getEdges()) {
                        if (!edge.isReference()) {
                            events.add(getUnresolvedEvent(edge.getInterval(), partitionIndex));
                        }
                    }
                }
            }

            /*
            final BiFunction<SVInterval, SVInterval, Boolean> repartitionFunction = (a, b) -> graphRepartitioningFunction(a, b, 0.1);
            final SVGraphPartitioner repartitioner = new SVGraphPartitioner(partition);
            final List<SVGraph> repartitions = repartitioner.getIndependentSubgraphs(repartitionFunction);
            for (final SVGraph repartition : repartitions) {
                final Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> retryResult = generateEvents(repartition, partitionIndex,
                        arguments.minEventProb, arguments.maxPathLengthFactor, arguments.maxEdgeVisits, copyNumberPosteriorsTree,
                        arguments.maxBranches, baselineCopyNumber, arguments.minEventSize, arguments.minHaplotypeProb, 1);
                if (retryResult == null) {
                    for (final IndexedSVGraphEdge edge : repartition.getEdges()) {
                        if (!edge.isReference()) {
                            events.add(getUnresolvedEvent(edge.getInterval(), partitionIndex));
                        }
                    }
                } else {
                    haplotypes.addAll(retryResult._1);
                    events.addAll(retryResult._2);
                    //TODO this breaks child baseline copy state
                    //partitionHaplotypes.addAll(retryResult._1);
                }
            }
            */
        } else {
            partitionHaplotypes.addAll(result._1);
            //TODO this list gets too big
            //haplotypes.addAll(result._1);
            events.addAll(result._2);
        }

        //TODO does not work with intra-chromosomal events
        final SVIntervalTree<IndexedSVGraphEdge> referenceEdgeTree = new SVIntervalTree<>();
        for (final IndexedSVGraphEdge referenceEdge : partition.getReferenceEdges()) {
            referenceEdgeTree.put(referenceEdge.getInterval(), referenceEdge);
        }
        for (final Integer childIndex : partitionDependenceGraph.get(partitionIndex)) {
            final SVGraph child = partitions.get(childIndex);
            final Double[] childCopyStateQualities = childCopyStates(child, partition.getReferenceEdges(), partitionHaplotypes, referenceEdgeTree, baselineCopyNumber);
            for (int copyState = 0; copyState < childCopyStateQualities.length; copyState++) {
                if (childCopyStateQualities[copyState] != null) {
                    final Tuple2<Collection<CalledSVGraphGenotype>, Collection<CalledSVGraphEvent>> childResult = callEventsHelper(childIndex, partitions, processedPartitions, copyState, childCopyStateQualities[copyState], partitionDependenceGraph);
                    processedPartitions[childIndex] = false; //Safe as long as there can be no graph loops
                    haplotypes.addAll(childResult._1);
                    events.addAll(childResult._2);
                }
            }
            processedPartitions[childIndex] = true;
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

    private static Double[] childCopyStates(final SVGraph child,
                                            final List<IndexedSVGraphEdge> referenceEdges,
                                            final List<CalledSVGraphGenotype> partitionHaplotypes,
                                            final SVIntervalTree<IndexedSVGraphEdge> referenceEdgeTree,
                                            final int baselineCopyNumber) {
        final int overlappingEdgeIndex = referenceEdgeTree.overlappers(child.getContigIntervals().get(0)).next().getValue().getIndex();
        /*
        int overlappingEdgeIndex = -1;
        //Assume child and graph are on the same contig
        final SVInterval childInterval = child.getContigIntervals().get(0); // TODO assumes not intrachromosomal
        final int pos = childInterval.getStart();
        for (int i = 0; i < referenceEdges.size(); i++) {
            final SVInterval edgeInterval = referenceEdges.get(i).getInterval();
            if (edgeInterval.getStart() <= pos && edgeInterval.getEnd() >= pos) {
                overlappingEdgeIndex = referenceEdges.get(i).getIndex();
                break;
            }
        }
        */

        final int MAX_COPY_STATE = 6;
        final Double[] childCopyStateQualities = new Double[MAX_COPY_STATE + 1];
        if (overlappingEdgeIndex != -1) {
            for (final CalledSVGraphGenotype genotype : partitionHaplotypes.stream().limit(100).collect(Collectors.toList())) { //TODO limiting genotypes
                int count = 0;
                for (final IndexedSVGraphPath haplotype : genotype.getHaplotypes()) {
                    for (final IndexedSVGraphEdge edge : haplotype.getEdges()) {
                        if (edge.getIndex() == overlappingEdgeIndex) {
                            count++;
                        }
                    }
                }
                if (count > MAX_COPY_STATE) {
                    throw new IllegalStateException("Genotype copy state was greater than the maximum of " + MAX_COPY_STATE);
                }
                if (childCopyStateQualities[count] == null || genotype.getProbability() < childCopyStateQualities[count]) {
                    childCopyStateQualities[count] = genotype.getProbability();
                }
            }
        } else {
            childCopyStateQualities[baselineCopyNumber] = Double.valueOf(Integer.MIN_VALUE);
        }
        return childCopyStateQualities;
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
        private final boolean homozygous;

        public SVGraphEvent(final CalledSVGraphEvent.Type type, final SVInterval interval,
                            final int groupId, final int pathId, final double probability,
                            final double evidenceProbability, final boolean resolved,
                            final boolean homozygous) {
            Utils.nonNull(type, "Type cannot be null");
            Utils.nonNull(interval, "Interval cannot be null");
            this.type = type;
            this.interval = interval;
            this.groupId = groupId;
            this.pathId = pathId;
            this.probability = probability;
            this.evidenceProbability = evidenceProbability;
            this.resolved = resolved;
            this.homozygous = homozygous;
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

        public boolean isHomozygous() { return homozygous; }
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
