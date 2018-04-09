package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import spire.std.SeqVectorEq;

import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.Stack;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Haplotype stream based on a {@link SeqGraph}.
 * <p>
 * We assume that there graph provided is a fully connected DAG with a single reference source and sink vertices
 * Where any path starting in the source ends up in the sink.
 * </p>
 *
 * The streamer provide a few way to stream through the haplotypes in the graph:
 * <ul>
 *     <li>{@link #bestScoreFirst()} exhaustive list of all haplotypes sorted by their {@code score}.</li>
 * </ul>
 */
public class SeqGraphHaplotypeStreamer {

    private final SeqGraph graph;
    private final SeqVertex source;
    private final SeqVertex sink;
    private final long numberOfHaplotypes;


    public SeqGraphHaplotypeStreamer(final SeqGraph graph) {
        Utils.nonNull(graph);
        this.graph = graph;
        this.source = graph.getReferenceSourceVertex();
        this.sink = graph.getReferenceSinkVertex();
        this.numberOfHaplotypes = calculateNumberOfHaplotypes(graph, source, sink);
    }

    private static long calculateNumberOfHaplotypes(final SeqGraph graph, final SeqVertex source, final SeqVertex sink) {
        final Map<SeqVertex, Long> results = new HashMap<>(graph.vertexSet().size());
        long total = 0;
        for (final SeqVertex v : graph.outgoingVerticesOf(source)) {
            total += calculateNumberOfHaplotypes(results, graph, v, sink);
        }
        return total;
    }

    private static long calculateNumberOfHaplotypes(final Map<SeqVertex, Long> results, SeqGraph graph, SeqVertex v, SeqVertex sink) {
        final Long knowResult = results.get(v);
        if (knowResult != null) {
            return knowResult;
        } else if (v == sink) {
            results.put(v, 1L);
            return 1L;
        } else {
            long result = 0;
            for (final SeqVertex vertex : graph.outgoingVerticesOf(v)) {
                result += calculateNumberOfHaplotypes(results, graph, vertex, sink);
            }
            results.put(v, result);
            return result;
        }
    }

    public long getNumberOfHaplotypes() {
        return numberOfHaplotypes;
    }

    public boolean veryLargeNumberOfHaplotypes() {
        return numberOfHaplotypes > 100;
    }

    /**
     * Returns a stream across all the possible haplotypes from the reference source
     * to sink vertices in a greedy deep-first fashion.
     * <p>
     *     Notice that in some instance for very complex graphs that may end up being a very long succession of haplotypes that
     *     don't cover varition across the graph well.
     * </p>
     * <p>
     *     However is a very efficient way to obtain an exhaustive list of haplotypes for relatively simple graphs.
     * </p>
     * @return never {@code null}
     */
    public Stream<Haplotype> deepFirst() {
        return StreamSupport.stream(new DeepFirstSpliterator(source), false);
    }

    /**
     * Returns a stream of haplotypes sorted by their scores.
     *
     * <p>
     *     Notice that in some instances some relevant variants might be missed amongst the top
     *     haplotypes if these are confined to sections of the graphs accessible thru a relative
     *     infrequent kmer transition (due to unequal allele balance furcation earlier in the graph).
     *     if the alternative path happens to have many possible paths.
     * </p>
     * @return never {@code null}
     */
    public Stream<Haplotype> bestScoreFirst() {
        final KBestHaplotypeFinder oldFinder = new KBestHaplotypeFinder(graph, source, sink);
        return Utils.stream(oldFinder).map(KBestHaplotype::haplotype)
                .map(h -> {
                    if (Double.isNaN(h.getScore())) {
                        h.setScore(oldFinder.score(h));
                    }
                    return h;
                });
    }

    /**
     * Generate random sample of haplotypes with replacement where the likelihood of each given
     * haplotype to be generated is proportional to the probability of its path from the
     * reference source vertex to the reference sink.
     *
     * <p>
     *     Notice that this Stream will return repeated haplotypes specially if there are just a few.
     * </p>
     *
     * @param rdn number generator to use for sampling.
     * @param distinct whether the returns sequence of haplotypes must be without replacement.
     * @param maxFailures in case distinct is true, the maximum number of failures to generate
     *    a new haplotype by likelihood sample before we give up.
     * @return never null.
     */
    public Stream<Haplotype> sampleByLikelihood(final Random rdn, final boolean distinct, final int maxFailures) {
        return StreamSupport.stream(new SampleByProbabilitySpliterator(rdn, distinct, maxFailures), true);
    }

    private Map<SeqVertex, Integer> calculateDensities() {
        return graph.vertexSet().stream().collect(Collectors.toMap(
                v -> v,
                v -> {
                    final Set<BaseEdge> incoming = graph.incomingEdgesOf(v);
                    final Set<BaseEdge> outgoing = graph.outgoingEdgesOf(v);
                    final int incomingSum = incoming.stream().mapToInt(BaseEdge::getMultiplicity).sum();
                    final int outgoingSum = outgoing.stream().mapToInt(BaseEdge::getMultiplicity).sum();
                    return Math.max(incomingSum, outgoingSum);
                }
        ));
    }

    private long[] calculateAccumulateDensity(final SeqVertex[] vertices, final Map<SeqVertex, Integer> densities) {
        final long[] result = new long[vertices.length];
        result[0] = densities.get(vertices[0]);
        for (int i = 1; i < vertices.length; i++) {
            result[i] = result[i - 1] + densities.get(vertices[i]);
        }
        return result;
    }

    private static class State {
        private final SeqVertex vertex;
        private final boolean isReference;
        private final double score;
        private final int prefixLength;

        private State(final SeqVertex vertex, final boolean isRef, final int prefixLength, final double score) {
            this.vertex = vertex;
            this.isReference = isRef;
            this.score = score;
            this.prefixLength = prefixLength;
        }
    }

    private class DeepFirstSpliterator implements Spliterator<Haplotype> {

        private long soFar = 0;
        private Stack<State> stack;
        final StringBuilder sequence;

        private DeepFirstSpliterator(final SeqVertex source) {
            super();
            sequence = new StringBuilder(1000);
            stack = new Stack<>();
            stack.push(new State(source, true, 0, 0));
        }

        @Override
        public boolean tryAdvance(Consumer<? super Haplotype> action) {
            while (!stack.isEmpty()) {
                final State next = stack.pop();
                sequence.setLength(next.prefixLength);
                sequence.append(next.vertex.getSequenceString());
                if (next.vertex == sink) {
                    final Haplotype haplotype;
                    try {
                        haplotype = new Haplotype(sequence.toString().getBytes(), next.isReference);
                    } catch (final RuntimeException ex) {
                        throw ex;
                    }
                    haplotype.setScore(next.score);
                    action.accept(haplotype);
                    soFar++;
                    return true;
                } else {
                    final Set<BaseEdge> edges = graph.outgoingEdgesOf(next.vertex);
                    int totalMultiplicity = 0;
                    for (final BaseEdge edge : edges) {
                        totalMultiplicity += edge.getMultiplicity();
                    }
                    final double totalMultiplicityDbl = totalMultiplicity;
                    for (final BaseEdge edge : edges) {
                        final SeqVertex child = graph.getEdgeTarget(edge);
                        final State childState =
                                new State(child, next.isReference && edge.isRef(), sequence.length(),
                                next.score + Math.log10(edge.getMultiplicity() / totalMultiplicityDbl));
                        stack.push(childState);
                    }
                }
            }
            return false;
        }

        @Override
        public Spliterator<Haplotype> trySplit() {
            return null;
        }

        @Override
        public long estimateSize() {
            return numberOfHaplotypes - soFar;
        }

        @Override
        public int characteristics() {
            return Spliterator.NONNULL & Spliterator.DISTINCT & Spliterator.SIZED;
        }
    }

    private class SampleByProbabilitySpliterator implements Spliterator<Haplotype> {
        private final Random rdn;
        private final boolean distinct;
        private final int maxFailures;
        private StringBuilder sequence;
        private final Set<Haplotype> soFar;

        public SampleByProbabilitySpliterator(final Random rdn, final boolean distinct, final int maxFailures) {
            this.rdn = Utils.nonNull(rdn);
            this.distinct = distinct;
            this.maxFailures = distinct ? ParamUtils.isPositiveOrZero(maxFailures, "the number of maximum failures must be 0 or greater") : maxFailures;
            this.sequence = new StringBuilder(1500);
            this.soFar = distinct ? new HashSet<>() : null;
        }

        private Haplotype nextHaplotype() {
            SeqVertex  current = source;
            sequence.setLength(0);
            double score = 0.0;
            boolean isRef = true;
            sequence.append(current.getSequenceString());
            while (current != sink) {
                final Set<BaseEdge> outgoingEdges = graph.outgoingEdgesOf(current);
                int totalWeight = 0;
                for (final BaseEdge edge : outgoingEdges) {
                    final int weight = Math.max(edge.getMultiplicity(), 1);
                    totalWeight += weight;
                }
                score -= Math.log10(totalWeight);
                final int accumulatedWeight = rdn.nextInt(totalWeight);
                for (final BaseEdge edge : outgoingEdges) {
                    final int weight = Math.max(edge.getMultiplicity(), 1);
                    totalWeight -= weight;
                    if (totalWeight <= accumulatedWeight) { // this is guaranteed to happen at some point.
                        current = graph.getEdgeTarget(edge);
                        sequence.append(current.getSequenceString());
                        isRef &= edge.isRef();
                        score += Math.log10(weight);
                        break;
                    }
                }
            }
            final Haplotype haplotype = new Haplotype(sequence.toString().getBytes(), isRef);
            haplotype.setScore(score);
            return haplotype;
        }

        @Override
        public boolean tryAdvance(final Consumer<? super Haplotype> action) {
            if (!distinct) {
                action.accept(nextHaplotype());
                return true;
            } else {
                int attempts = 0;
                do {
                  final Haplotype haplotype = nextHaplotype();
                  if (soFar.add(haplotype)) {
                      action.accept(haplotype);
                      return true;
                  }
                  attempts++;
                } while (attempts < maxFailures);
                return false;
            }
        }

        @Override
        public Spliterator<Haplotype> trySplit() {
            if (!distinct) {
                return new SampleByProbabilitySpliterator(new Random(rdn.nextInt()), false, 0);
            } else {
                return null; // not support for parelism if asked for distinct output.
            }
        }

        @Override
        public long estimateSize() {
            return Long.MAX_VALUE;
        }

        @Override
        public int characteristics() {
            return Spliterator.NONNULL;
        }
    }
}
