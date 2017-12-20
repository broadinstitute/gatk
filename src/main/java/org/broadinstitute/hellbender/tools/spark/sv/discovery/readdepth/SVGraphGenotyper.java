package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Enumerates genotypes on an SV graph
 */
public final class SVGraphGenotyper {

    final SVGraph graph;

    public SVGraphGenotyper(final SVGraph graph) {
        this.graph = graph;
    }

    /**
     * Normalizes genotype likelihoods
     */
    private static final void setGenotypedProbabilities(final Collection<SVGraphGenotype> genotypes) {
        final double logDenom = Math.log(genotypes.stream().mapToDouble(SVGraphGenotype::getLikelihood)
                .map(Math::exp)
                .map(val -> Math.max(val, Double.MIN_NORMAL))
                .sum());
        for (final SVGraphGenotype genotype : genotypes) {
            genotype.setProbability(Math.exp(genotype.getLikelihood() - logDenom));
        }
    }

    /**
     * Edge visit for graph traversal
     */
    private static final void visitEdge(final IndexedSVGraphEdge edge, final SVGraphPathState path, final Queue<SVGraphPathState> queue, final int nodeIndex, final int maxEdgeVisits) {
        final IndexedSVGraphEdge augmentingEdge;
        if (edge.causesStrandSwitch()) {
            augmentingEdge = path.isCurrentlyInverted() ? edge.invertedCopy().copyWithUpstreamNodeAs(nodeIndex) : edge.copyWithUpstreamNodeAs(nodeIndex);
        } else {
            augmentingEdge = path.isCurrentlyInverted() ? edge.invertedCopy() : edge.nonInvertedCopy();
        }
        final SVGraphPathState augmentedPath = path.copy();
        if (augmentedPath.addEdge(augmentingEdge, maxEdgeVisits)) {
            if (augmentingEdge.causesStrandSwitch()) {
                augmentedPath.invert();
            }
            queue.add(augmentedPath);
        }
    }

    /**
     * Creates/adds haplotype when the path has reached an end node
     */
    private static final void completeHaplotype(final SVGraphPathState path, final Set<SVGraphGenotype> genotypes, final Queue<SVGraphPathState> queue, final Set<Integer> startingNodes, final int groupId, final int baselineCopyNumber) {
        final int currentPloidy = path.getCurrentPloidy();
        if (currentPloidy == baselineCopyNumber) {
            //Finished all haplotypes
            final double copyNumberScore = path.logLikelihood();
            final List<IndexedSVGraphPath> haplotypePaths = path.getPathsList().stream().map(edgesList -> new IndexedSVGraphPath(edgesList)).collect(Collectors.toList());
            final SVGraphGenotype genotype = new SVGraphGenotype(groupId, genotypes.size(), copyNumberScore, haplotypePaths);
            genotypes.add(genotype);
        } else if (currentPloidy < baselineCopyNumber) {
            //Begin new haplotype
            for (final Integer startingNode : startingNodes) {
                final SVGraphPathState augmentedPath = path.copy();
                augmentedPath.setCurrentNode(startingNode);
                augmentedPath.incrementPloidy();
                queue.add(augmentedPath);
            }
        }
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

    private List<SVCopyNumberInterval> getSortedOverlappingCopyNumberIntervals(final Collection<IndexedSVGraphEdge> edges, final SVIntervalTree<SVCopyNumberInterval> copyNumberPosteriorsTree) {
        return edges.stream()
                .flatMap(edge -> Utils.stream(copyNumberPosteriorsTree.overlappers(edge.getInterval())))
                .map(SVIntervalTree.Entry::getValue)
                .distinct()
                .sorted(SVIntervalUtils.getCopyNumberIntervalDictionaryOrderComparator())
                .collect(Collectors.toList());
    }

    /**
     * Gets copy number posteriors for all edges
     */
    private List<EdgeCopyNumberPosterior> getEdgeCopyNumberPosteriors(final SVIntervalTree<SVCopyNumberInterval> copyNumberPosteriorsTree) {

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
     * Enumerates genotypes with breadth-first search
     */
    public Collection<SVGraphGenotype> enumerate(final SVIntervalTree<SVCopyNumberInterval> copyNumberPosteriorsTree,
                                                 final int groupId, final int baselineCopyNumber,
                                                 final double maxPathLengthFactor, final int maxEdgeVisits,
                                                 final int maxQueueSize) {

        final List<SVGraphNode> nodes = graph.generateNodes();
        final Set<Integer> startingNodes = graph.getStartingNodes();
        final Set<Integer> endingNodes = graph.getEndingNodes();

        final List<EdgeCopyNumberPosterior> copyNumberPosteriors = getEdgeCopyNumberPosteriors(copyNumberPosteriorsTree);
        final int maxPathLength = (int) (maxPathLengthFactor * graph.getReferenceEdges().size());

        final Set<SVGraphGenotype> genotypes = new HashSet<>();
        final Queue<SVGraphPathState> queue = new ArrayDeque<>();
        for (final Integer startingNode : startingNodes) {
            queue.add(new SVGraphPathState(startingNode, copyNumberPosteriors));
        }
        while (!queue.isEmpty()) {
            if (queue.size() > maxQueueSize) return null;
            final SVGraphPathState path = queue.poll();
            final Integer nodeIndex = path.getCurrentNode();

            //Cap max path length (extremely long paths caused by degenerate loops)
            if (path.getCurrentPath().size() > maxPathLength) continue;

            //Finish haplotype if we reach the other end of the event in an uninverted state
            if (!path.isCurrentlyInverted() && endingNodes.contains(nodeIndex)) {
                completeHaplotype(path, genotypes, queue, startingNodes, groupId, baselineCopyNumber);
            }

            //Visit neighbors
            final boolean lastEdgeWasReference = path.getCurrentPath().isEmpty() || path.getCurrentPath().get(path.getCurrentPath().size() - 1).isReference();
            final SVGraphNode node = nodes.get(nodeIndex);
            final List<IndexedSVGraphEdge> neighborEdges = path.isCurrentlyInverted() ? node.getInEdges() : node.getOutEdges();
            for (final IndexedSVGraphEdge edge : neighborEdges) {
                //Breakpoint double-jumps are invalid
                if (edge.isReference() || lastEdgeWasReference) {
                    visitEdge(edge, path, queue, nodeIndex, maxEdgeVisits);
                }
            }
        }
        setGenotypedProbabilities(genotypes);
        return genotypes;
    }

    /**
     * Container for the state of a genotype during graph traversal
     */
    private static final class SVGraphPathState {
        private final List<List<IndexedSVGraphEdge>> pathsList;
        private final List<EdgeCopyNumberPosterior> edgeCopyNumberPosteriors;
        private final int[] edgeCopyNumberStates;
        private final int numCopyNumberStates;
        private int currentNode;
        private int currentPloidy;
        boolean currentlyInverted;

        public SVGraphPathState(final int initialNode, final List<EdgeCopyNumberPosterior> edgeCopyNumberPosteriors) {
            this.pathsList = new ArrayList<>();
            this.pathsList.add(new ArrayList<>());
            this.edgeCopyNumberPosteriors = edgeCopyNumberPosteriors;
            this.edgeCopyNumberStates = new int[edgeCopyNumberPosteriors.size()];
            this.numCopyNumberStates = edgeCopyNumberPosteriors.isEmpty() ? 0 : edgeCopyNumberPosteriors.get(0).numCopyNumberStates();
            this.currentNode = initialNode;
            currentPloidy = 1;
            currentlyInverted = false;
        }

        private SVGraphPathState(final int currentNode, final List<EdgeCopyNumberPosterior> edgeCopyNumberPosteriors,
                                 final int[] edgeCopyNumberStates, final int currentPloidy,
                                 final boolean currentlyInverted, final int numCopyNumberStates) {
            this.currentNode = currentNode;
            this.pathsList = new ArrayList<>();
            this.edgeCopyNumberPosteriors = edgeCopyNumberPosteriors;
            this.edgeCopyNumberStates = Arrays.copyOf(edgeCopyNumberStates, edgeCopyNumberStates.length);
            this.currentPloidy = currentPloidy;
            this.currentlyInverted = currentlyInverted;
            this.numCopyNumberStates = numCopyNumberStates;
        }

        public boolean isCurrentlyInverted() {
            return currentlyInverted;
        }

        public void invert() {
            currentlyInverted = !currentlyInverted;
        }

        public List<List<IndexedSVGraphEdge>> getPathsList() {
            return pathsList;
        }

        public List<IndexedSVGraphEdge> getCurrentPath() {
            return pathsList.get(currentPloidy - 1);
        }

        public int getCurrentNode() {
            return currentNode;
        }

        public void setCurrentNode(final int currentNode) {
            this.currentNode = currentNode;
        }

        public double logLikelihood() {
            return IntStream.range(0, edgeCopyNumberPosteriors.size()).mapToDouble(k -> edgeCopyNumberPosteriors.get(k).getLogPosterior(edgeCopyNumberStates[k])).sum();
        }

        public boolean addEdge(final IndexedSVGraphEdge edge, final double maxState) {
            getCurrentPath().add(edge);
            if (edge.isReference()) {
                edgeCopyNumberStates[edge.getIndex()] += 1;
                if (edgeCopyNumberStates[edge.getIndex()] > maxState) return false;
            }
            currentNode = edge.getDownstreamNodeIndex();
            if (currentNode == -1) {
                throw new GATKException("Attempted to add edge without a set upstream node");
            }
            return true;
        }

        public void incrementPloidy() {
            currentPloidy++;
            pathsList.add(new ArrayList<>());
        }

        public int getCurrentPloidy() {
            return currentPloidy;
        }

        public SVGraphPathState copy() {
            final SVGraphPathState copy = new SVGraphPathState(currentNode, edgeCopyNumberPosteriors, edgeCopyNumberStates, currentPloidy, currentlyInverted, numCopyNumberStates);
            for (final List<IndexedSVGraphEdge> path : pathsList) {
                final List<IndexedSVGraphEdge> pathCopy = new ArrayList<>(path.size());
                for (final IndexedSVGraphEdge edge : path) {
                    pathCopy.add(edge);
                }
                copy.pathsList.add(pathCopy);
            }
            return copy;
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
