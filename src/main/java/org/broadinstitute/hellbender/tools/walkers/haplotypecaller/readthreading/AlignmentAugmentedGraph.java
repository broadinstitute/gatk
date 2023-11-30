package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.api.client.util.Lists;
import com.google.api.client.util.Sets;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang.math.IntRange;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.ChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class AlignmentAugmentedGraph {

    // TODO: magic constant, make adjustable
    private static final int MAX_BRANCH_VERTICES = 100;

    private final int kmerSize;

    private final byte minBaseQualityToUseInAssembly;

    public AlignmentAugmentedGraph(final int kmerSize, final byte minBaseQualityToUseInAssembly) {
        this.kmerSize = kmerSize;
        this.minBaseQualityToUseInAssembly = minBaseQualityToUseInAssembly;
    }

    public void doAssembly(final Haplotype refHaplotype, final Iterable<GATKRead> reads, final ChainPruner pruner) {
        // make an prune a simple de Bruijn graph in order to know which kmers are worth tracking
        final PlainDeBruijnGraph initialGraph = new PlainDeBruijnGraph(kmerSize, minBaseQualityToUseInAssembly);
        initialGraph.addSequence("ref", refHaplotype.getBases(), 1, true);
        reads.forEach(read -> initialGraph.addRead(read, null));
        initialGraph.buildGraphIfNecessary();
        pruner.pruneLowWeightChains(initialGraph);

        final Set<Kmer> kmersAfterPruning = initialGraph.vertexSet().stream()
                .map(vertex -> new Kmer(vertex.getSequence()))
                .collect(Collectors.toSet());

        final List<Pair<Kmer, IntRange>> rangedKmers = Utils.stream(reads).flatMap(read -> kmerizeRead(read).stream())
                .filter(pair -> pair.getLeft() != null).toList();

        final List<Pair<Kmer, Integer>> unambiguousKmers = rangedKmers.stream()
                .filter(pair -> kmersAfterPruning.contains(pair.getLeft()))
                .filter(pair -> Math.abs(pair.getRight().getMaximumInteger() - pair.getRight().getMinimumInteger()) < kmerSize)
                .map(pair -> Pair.of(pair.getLeft(), (pair.getRight().getMinimumInteger() + pair.getRight().getMaximumInteger())/2))
                .toList();

        final VertexManager vertexManager = new VertexManager(kmerSize, unambiguousKmers);

        final AugmentedKmerGraph graph = makeAugmentedKmerGraph(reads, vertexManager);

        final List<AugmentedVertex> branchVertices = graph.vertexSet().stream()
                .filter(v -> graph.outDegreeOf(v) > 1).toList();

        // decision vertices = anything after a branch or a source (wince which source to begin with is a decision)
        final List<AugmentedVertex> decisionVertices = Stream.concat(graph.getSources().stream(),
                branchVertices.stream().flatMap(bv -> graph.outgoingVerticesOf(bv).stream()))
                .distinct() // needed because a decision vertex may have multiplo parents
                .toList();

        final Map<AugmentedVertex, Integer> decisionVertexIndexMap = IntStream.range(0, decisionVertices.size()).boxed()
                .collect(Collectors.toMap(decisionVertices::get, n->n));

        final List<Integer> sourceVertexIndices = graph.getSources().stream()
                .map(decisionVertexIndexMap::get).toList();
        
        final Map<GATKRead, int[]> readFeatureMap = Utils.stream(reads)
                .collect(Collectors.toMap(read -> read, read -> featurizeRead(vertexManager, graph, decisionVertices, decisionVertexIndexMap, sourceVertexIndices, read)));

    }

    private int[] featurizeRead(final VertexManager vertexManager, final AugmentedKmerGraph graph, final List<AugmentedVertex> decisionVertices,
                                final Map<AugmentedVertex, Integer> decisionVertexIndexMap, final List<Integer> sourceVertexIndices, GATKRead read) {
        AugmentedVertex lastVertexInRead = null;
        final int[] features = new int[decisionVertices.size()];

        for (final Pair<Kmer, IntRange> pair : kmerizeRead(read)) {
            final Kmer kmer = pair.getLeft();
            final IntRange positionRange = pair.getRight();
            final Optional<AugmentedVertex> maybeVertex = kmer == null ? Optional.empty() : vertexManager.getVertex(kmer, positionRange);

            if (maybeVertex.isPresent()) {
                final AugmentedVertex vertex = maybeVertex.get();
                // if last vertex was null due to being at the read start or bad bases disqualifying a kmer in
                // the middle of the read, traverse backwards until the last decision vertex
                if (lastVertexInRead == null) {
                    AugmentedVertex backwardVertex = vertex;
                    while (!decisionVertexIndexMap.containsKey(backwardVertex)) {
                        backwardVertex = graph.getEdgeSource(graph.incomingEdgeOf(backwardVertex));
                    }

                    if (graph.isSource(backwardVertex)) {
                        sourceVertexIndices.forEach(n -> features[n] = -1);
                        features[decisionVertexIndexMap.get(backwardVertex)] = 1;
                    } else if (graph.inDegreeOf(backwardVertex) == 1) { // in rare edge case of two parents, it's ambiguous and we leave the feature at zero
                        final AugmentedVertex parent = graph.getEdgeSource(graph.incomingEdgeOf(backwardVertex));
                        graph.outgoingVerticesOf(parent).stream()
                                .filter(decisionVertexIndexMap::containsKey)
                                .map(decisionVertexIndexMap::get)
                                .forEach(n -> features[n] = -1);
                        features[decisionVertexIndexMap.get(backwardVertex)] = 1;
                    }
                } else if (decisionVertexIndexMap.containsKey(vertex)) {  // decision vertex
                    // set sibling features to -1 and this feature to +1
                    graph.outgoingVerticesOf(lastVertexInRead).forEach(v -> features[decisionVertexIndexMap.get(v)] = -1);
                    features[decisionVertexIndexMap.get(vertex)] = 1;
                }
            }

            lastVertexInRead = maybeVertex.orElse(null);
        }
        return features;
    }

    private AugmentedKmerGraph makeAugmentedKmerGraph(Iterable<GATKRead> reads, VertexManager vertexManager) {
        final AugmentedKmerGraph graph = new AugmentedKmerGraph(kmerSize);
        vertexManager.allVertices().forEach(graph::addVertex);

        // TODO: thread the reference sequence???
        for (final GATKRead read : reads) {
            final List<Pair<Kmer, IntRange>> kmers = kmerizeRead(read);
            for (int n = 0; n < kmers.size() - 1; n++) {
                final Kmer kmer1 = kmers.get(n).getLeft();
                final Kmer kmer2 = kmers.get(n+1).getLeft();

                if (kmer1 != null && kmer2 != null) {   // null if unusable kmers due to 'N' or low-qual bases
                    final Optional<AugmentedVertex> vertex1 = vertexManager.getVertex(kmer1, kmers.get(n).getRight());
                    final Optional<AugmentedVertex> vertex2 = vertexManager.getVertex(kmer1, kmers.get(n+1).getRight());

                    if (vertex1.isPresent() && vertex2.isPresent()) {
                        graph.addEdge(vertex1.get(), vertex2.get());
                    }
                }
            }
        }

        graph.removeSingletonOrphanVertices();
        return graph;
    }

    // extracts kmers and their alignment positions -- soft clips are assigned a range with one ambiguous end eg
    // [-infinity, 10] or [10, infinity]
    private List<Pair<Kmer, IntRange>> kmerizeRead(final GATKRead read) {
        final byte[] sequence = read.getBases();
        final byte[] qualities = read.getBaseQualities();
        final List<Pair<Kmer, IntRange>> result = new ArrayList<>(sequence.length - kmerSize + 1);

        final int start = read.hasAttribute(ReadUtils.ORIGINAL_SOFTCLIP_START_TAG) ?
                read.getAttributeAsInteger(ReadUtils.ORIGINAL_SOFTCLIP_START_TAG) : read.getStart();
        final int end = read.hasAttribute(ReadUtils.ORIGINAL_SOFTCLIP_END_TAG) ?
                read.getAttributeAsInteger(ReadUtils.ORIGINAL_SOFTCLIP_END_TAG) : read.getEnd();

        int lastUnusableReadOffset = -1;   // last kmer-disqualifying base within the current kmer span
        for (int n = Math.min(kmerSize, sequence.length) - 1; n >= 0; n--) {
            if (!baseIsUsableForAssembly(sequence[n], qualities[n])) {
                lastUnusableReadOffset = n;
                break;
            }
        }

        int readOffset = 0;
        int refPosition = read.getSoftStart();
        for (final CigarElement cigarElement : read.getCigarElements()) {
            final CigarOperator op = cigarElement.getOperator();
            final int length = cigarElement.getLength();

            if (!op.consumesReadBases()) {  // deletion: skip ahead on reference but don't make kmers
                refPosition += op.consumesReferenceBases() ? length : 0;
            } else {    // make kmers
                for (int n = 0; n < length; n++) {
                    if (readOffset + kmerSize > sequence.length) {  // kmer would go off end of read
                        break;
                    } else {
                        final int endOfKmer = readOffset + kmerSize - 1;
                        if (!baseIsUsableForAssembly(sequence[endOfKmer], qualities[endOfKmer])) {
                            lastUnusableReadOffset = endOfKmer;
                        }

                        if (lastUnusableReadOffset < readOffset) {  // make a kmer
                            final Kmer kmer = new Kmer(sequence, readOffset, kmerSize);
                            final int rangeMin = refPosition < start ? Integer.MIN_VALUE : refPosition;
                            // TODO: verify that this is > and not >=
                            final int rangeMax = refPosition > end ? Integer.MAX_VALUE : refPosition;
                            result.add(Pair.of(kmer, new IntRange(rangeMin, rangeMax)));
                        } else {
                            result.add(Pair.of(null, new IntRange(refPosition)));   // unusable kmer
                        }

                        refPosition += op.consumesReferenceBases() ? 1 : 0;
                        readOffset++;
                    }
                }
            }
        }
        return result;
    }

    private static class VertexManager {
        private final Map<Kmer, PositionClusterer> kmerMap = new HashMap<>();

        public VertexManager(final int tolerance, final Iterable<Pair<Kmer, Integer>> kmersAndPositions) {
            for (final Pair<Kmer, Integer> pair : kmersAndPositions) {
                kmerMap.computeIfAbsent(pair.getLeft(), kmer -> new PositionClusterer(tolerance, kmer)).add(pair.getRight());
            }
            kmerMap.values().forEach(clusterer -> clusterer.finish());
        }

        public Optional<AugmentedVertex> getVertex(final Kmer kmer, final IntRange range) {
            final int rangeMin = range.getMinimumInteger();
            final int rangeMax = range.getMaximumInteger();

            // TODO: magic conventions!
            if (rangeMin < 0) {
                return getVertexFromMaxPosition(kmer, rangeMax);
            } else if (rangeMax == Integer.MAX_VALUE) {
                return getVertexFromMinPosition(kmer, rangeMin);
            } else {
                return getVertex(kmer, (rangeMin + rangeMax)/2);
            }
        }

        public Optional<AugmentedVertex> getVertex(final Kmer kmer, final int position) {
            final PositionClusterer clusterer = kmerMap.get(kmer);
            return clusterer == null ? Optional.empty() : clusterer.getVertex(position);
        }

        public Optional<AugmentedVertex> getVertexFromMaxPosition(final Kmer kmer, final int maxPosition) {
            final PositionClusterer clusterer = kmerMap.get(kmer);
            return clusterer == null ? Optional.empty() : clusterer.getVertexFromMaxPosition(maxPosition);
        }

        public Optional<AugmentedVertex> getVertexFromMinPosition(final Kmer kmer, final int minPosition) {
            final PositionClusterer clusterer = kmerMap.get(kmer);
            return clusterer == null ? Optional.empty() : clusterer.getVertexFromMinPosition(minPosition);
        }

        public Stream<AugmentedVertex> allVertices() {
            return kmerMap.values().stream().flatMap(clusterer -> clusterer.getVertices().stream());
        }

    }

    private static class PositionClusterer {
        private final Kmer kmer;

        // this will always get sorted from greatest to least count whenever thing get too wrong
        private final List<Pair<Integer, MutableInt>> clustersAndCounts = new ArrayList<>(1);

        private final int tolerance;

        private boolean finalized;

        private final List<AugmentedVertex> vertices = new ArrayList<>(1);


        public PositionClusterer(final int tolerance, final Kmer kmer) {
            this.tolerance = tolerance;
            this.kmer = kmer;
        }

        public List<AugmentedVertex> getVertices() { return vertices; }

        public void finish() {
            if (finalized) {
                return;
            }
            finalized = true;
            sort();
            clustersAndCounts.forEach(pair -> vertices.add(new AugmentedVertex(kmer.bases(), pair.getLeft())));
        }

        public Optional<AugmentedVertex> getVertex(final int position) {
            if (vertices.size() == 1) {
                return Optional.of(vertices.get(0));
            } else {
                for (int n = 0; n < clustersAndCounts.size(); n++) {
                    if (inCluster(position, clustersAndCounts.get(n).getLeft())) {
                        return Optional.of(vertices.get(n));
                    }
                }
            }
            return Optional.empty();
        }

        public Optional<AugmentedVertex> getVertexFromMaxPosition(final int maxPosition) {
            int matchIndex = -1;
            for (int n = 0; n < clustersAndCounts.size(); n++) {
                if (clustersAndCounts.get(n).getLeft() >= maxPosition) {
                    if (matchIndex != -1) {
                        return Optional.empty();    // multiple matches, ambiguous soft clip
                    }
                    matchIndex = n;
                }
            }

            return matchIndex != -1 ? Optional.of(vertices.get(matchIndex)) : Optional.empty();
        }

        public Optional<AugmentedVertex> getVertexFromMinPosition(final int minPosition) {
            int matchIndex = -1;
            for (int n = 0; n < clustersAndCounts.size(); n++) {
                if (clustersAndCounts.get(n).getLeft() >= minPosition) {
                    if (matchIndex != -1) {
                        return Optional.empty();    // multiple matches, ambiguous soft clip
                    }
                    matchIndex = n;
                }
            }

            return matchIndex != -1 ? Optional.of(vertices.get(matchIndex)) : Optional.empty();
        }



        public void add(final int position) {
            Utils.validate(!finalized, "already finalized");
            boolean newCluster = true;
            for (int n = 0; n < clustersAndCounts.size(); n++) {
                final MutableInt count = clustersAndCounts.get(n).getRight();
                if (inCluster(position, clustersAndCounts.get(n).getLeft())) {
                    count.increment();
                    if (n > 0 && count.getValue() > 2 * clustersAndCounts.get(0).getRight().getValue()) {
                        sort();
                    }
                    newCluster = false;
                    break;
                }
            }

            if (newCluster) {
                clustersAndCounts.add(Pair.of(position, new MutableInt(1)));
            }
        }

        private void sort() {
            clustersAndCounts.sort(Comparator.comparingInt(pair -> -pair.getRight().getValue()));
        }

        private boolean inCluster(final int position, final int cluster) {
            return Math.abs(position - cluster) < tolerance;
        }
    }

    protected boolean baseIsUsableForAssembly(final byte base, final byte qual) {
        return base != BaseUtils.Base.N.base && qual >= minBaseQualityToUseInAssembly;
    }
}
