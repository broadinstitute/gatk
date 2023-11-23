package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.api.client.util.Lists;
import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang.math.IntRange;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.ChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class AlignmentAugmentedGraph extends AbstractReadThreadingGraph {

    // TODO: magic constant, make adjustable
    private static final int MAX_BRANCH_VERTICES = 100;

    // TODO: magic constant, make adjustable
    private static final int MAX_CLUSTERS_PER_KMER = 5;

    public AlignmentAugmentedGraph(final int kmerSize, final boolean debugGraphTransformations, final byte minBaseQualityToUseInAssembly) {
        super(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly, 1, -1);
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

        final AugmentedKmerGraph graph = new AugmentedKmerGraph(kmerSize);
        vertexManager.allVertices.forEach(graph::addVertex);

        // TODO: thread the reference sequence???
        for (final GATKRead read : reads) {
            final List<Pair<Kmer, IntRange>> kmers = kmerizeRead(read);
            for (int n = 0; n < kmers.size() - 1; n++) {
                final Kmer kmer1 = kmers.get(n).getLeft();
                final Kmer kmer2 = kmers.get(n+1).getLeft();

                if (kmer1 != null && kmer2 != null) {   // null if unusable kmers due to 'N' or low-qual bases
                    final Optional<MultiDeBruijnVertex> vertex1 = vertexManager.getVertex(kmer1, kmers.get(n).getRight());
                    final Optional<MultiDeBruijnVertex> vertex2 = vertexManager.getVertex(kmer1, kmers.get(n+1).getRight());

                    if (vertex1.isPresent() && vertex2.isPresent()) {
                        graph.addEdge(vertex1.get(), vertex2.get());
                    }
                }
            }
        }

        graph.removeSingletonOrphanVertices();




        final Set<Kmer> decisionKmers = findDecisionKmers(refHaplotype, reads, pruner);
        final Map<Kmer, TreeMultiset<Integer>> decisionKmerPositions = gatherDecisionKmerPositions(decisionKmers, rangedKmers);
        final Map<Kmer, List<Integer>> decisionKmerClusters = clusterDecisionKmers(decisionKmers, decisionKmerPositions);


    }

    private Set<Kmer> findDecisionKmers(final Haplotype refHaplotype, final Iterable<GATKRead> reads, final ChainPruner pruner) {
        // make a regular read threading graph in order to identify decision vertices
        final PlainDeBruijnGraph initialGraph = new PlainDeBruijnGraph(kmerSize, minBaseQualityToUseInAssembly);
        initialGraph.addSequence("ref", refHaplotype.getBases(), 1, true);
        reads.forEach(read -> initialGraph.addRead(read, null));
        initialGraph.buildGraphIfNecessary();
        pruner.pruneLowWeightChains(initialGraph);

        final List<MultiDeBruijnVertex> branchVertices = initialGraph.vertexSet().stream()
                .filter(v -> initialGraph.outDegreeOf(v) > 1).toList();

        final List<Pair<MultiDeBruijnVertex, Integer>> decisionVertices = new ArrayList<>();
        for (final MultiDeBruijnVertex vertex : branchVertices) {
            final List<MultiSampleEdge> edges = Lists.newArrayList(initialGraph.outgoingEdgesOf(vertex));
            final int[] multiplicities = edges.stream().mapToInt(MultiSampleEdge::getMultiplicity).toArray();
            final int branchiness = (int) MathUtils.sum(multiplicities) - MathUtils.arrayMax(multiplicities);
            // TODO: if there are more than 2 outgoing edges this gives minor edges more branchiness than they deserve
            edges.forEach(edge -> decisionVertices.add(Pair.of(initialGraph.getEdgeTarget(edge), branchiness)));
        }

        return decisionVertices.stream()
                .sorted(Comparator.comparingInt(Pair<MultiDeBruijnVertex, Integer>::getRight).reversed())
                .limit(MAX_BRANCH_VERTICES)
                .map(Pair::getLeft)
                .map(vertex -> new Kmer(vertex.getSequence()))
                .collect(Collectors.toSet());
    }

    private Map<Kmer, TreeMultiset<Integer>> gatherDecisionKmerPositions(final Set<Kmer> decisionKmers, final List<Pair<Kmer, IntRange>> rangedKmers) {
        final Map<Kmer, TreeMultiset<Integer>> decisionKmerPositions = new HashMap<>();
        decisionKmers.forEach(kmer -> decisionKmerPositions.put(kmer, TreeMultiset.create()));
        rangedKmers.stream()
                .filter(pair -> decisionKmers.contains(pair.getLeft())) // it's a decision kmer
                .filter(pair -> pair.getRight().getMaximumInteger() - pair.getRight().getMinimumInteger() < kmerSize)   // its location is roughly determined
                .forEach(pair -> decisionKmerPositions.get(pair.getLeft()).add(pair.getRight().getMinimumInteger()));
        return decisionKmerPositions;
    }

    private Map<Kmer, List<Integer>> clusterDecisionKmers(Set<Kmer> decisionKmers, Map<Kmer, TreeMultiset<Integer>> decisionKmerPositions) {
        final Map<Kmer, List<Integer>> decisionKmerClusters = new HashMap<>();
        decisionKmers.forEach(kmer -> decisionKmerClusters.put(kmer, Lists.newArrayList()));
        for (final Kmer kmer : decisionKmers) {
            final List<Pair<Integer, Integer>> clusters = Lists.newArrayList(); // positions and multiplicity
            final TreeMultiset<Integer> positions = decisionKmerPositions.get(kmer);
            int currentCluster = positions.firstEntry().getElement();
            int currentMultiplicity = 0;
            for (final Multiset.Entry<Integer> entry : positions.entrySet()) { // unique values is in ascending order
                // TODO: is this tolerance the right criterion?
                if (entry.getElement() < currentCluster + kmerSize) {   // expand current cluster if within tolerance
                    currentMultiplicity += positions.count(entry.getCount());
                } else {
                    clusters.add(Pair.of(currentCluster, currentMultiplicity));
                    currentCluster = entry.getElement();
                    currentMultiplicity = entry.getCount();
                }
            }
            clusters.add(Pair.of(currentCluster, currentMultiplicity)); // flush the last cluster
            final List<Integer> clusterList = decisionKmerClusters.get(kmer);
            clusters.stream().sorted(Comparator.comparingInt(pair -> -pair.getRight()))  // sort from greatest to least multiplciity
                    .map(Pair::getLeft).limit(MAX_CLUSTERS_PER_KMER).forEach(clusterList::add);
        }
        return decisionKmerClusters;
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

    /**
     * Checks whether a kmer can be the threading start based on the current threading start location policy.
     *
     * @param kmer the query kmer.
     * @return {@code true} if we can start thread the sequence at this kmer, {@code false} otherwise.
     * @see #setThreadingStartOnlyAtExistingVertex(boolean)
     */
    protected abstract boolean isThreadingStart(final Kmer kmer, final boolean startThreadingOnlyAtExistingVertex);

    // get the next kmerVertex for ChainExtension and validate if necessary.
    protected abstract MultiDeBruijnVertex getNextKmerVertexForChainExtension(final Kmer kmer, final boolean isRef, final MultiDeBruijnVertex prevVertex);

    // perform any necessary preprocessing on the graph (such as non-unique kmer determination) before the graph is constructed
    protected void preprocessReads() {
        throw new UnsupportedOperationException("HAven't implemented yet.");
    }

    // heuristic to decide if the graph should be thrown away based on low complexity/poor assembly
    @Override
    public boolean isLowQualityGraph() { return false; }

    // whether reads are needed after graph construction
    @Override
    protected boolean shouldRemoveReadsAfterGraphConstruction() {
        return false;
    }

    // Method that will be called immediately before haplotype finding in the event there are alteations that must be made to the graph based on implementation
    @Override
    public void postProcessForHaplotypeFinding(File debugGraphOutputPath, Locatable refHaplotype) {
        throw new UnsupportedOperationException("HAven't implemented yet.");
    }

    /**
     * Define the behavior for how the graph should keep track of a potentially new kmer.
     *
     * @param kmer      (potentially) new kmer to track
     * @param newVertex corresponding vertex for that kmer
     */
    protected abstract void trackKmer(Kmer kmer, MultiDeBruijnVertex newVertex);

    @Override
    public void recoverDanglingTails(final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll,
                                     final SmithWatermanAligner aligner, final SWParameters danglingTailSWParameters) {
        throw new UnsupportedOperationException("Not implemented yet but it's going to be different.");
    }

    @Override
    public void recoverDanglingHeads(final int pruneFactor, final int minDanglingBranchLength, final boolean recoverAll,
                                     final SmithWatermanAligner aligner, final SWParameters danglingTailSWParameters) {
        throw new UnsupportedOperationException("Not implemented yet but it's going to be different.");
    }

    private static class VertexManager {
        private final Map<Kmer, PositionClusterer> kmerMap = new HashMap<>();

        public VertexManager(final int tolerance, final Iterable<Pair<Kmer, Integer>> kmersAndPositions) {
            for (final Pair<Kmer, Integer> pair : kmersAndPositions) {
                kmerMap.computeIfAbsent(pair.getLeft(), kmer -> new PositionClusterer(tolerance, kmer)).add(pair.getRight());
            }
            kmerMap.values().forEach(clusterer -> clusterer.finish());
        }

        public Optional<MultiDeBruijnVertex> getVertex(final Kmer kmer, final IntRange range) {
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

        public Optional<MultiDeBruijnVertex> getVertex(final Kmer kmer, final int position) {
            final PositionClusterer clusterer = kmerMap.get(kmer);
            return clusterer == null ? Optional.empty() : clusterer.getVertex(position);
        }

        public Optional<MultiDeBruijnVertex> getVertexFromMaxPosition(final Kmer kmer, final int maxPosition) {
            final PositionClusterer clusterer = kmerMap.get(kmer);
            return clusterer == null ? Optional.empty() : clusterer.getVertexFromMaxPosition(maxPosition);
        }

        public Optional<MultiDeBruijnVertex> getVertexFromMinPosition(final Kmer kmer, final int minPosition) {
            final PositionClusterer clusterer = kmerMap.get(kmer);
            return clusterer == null ? Optional.empty() : clusterer.getVertexFromMinPosition(minPosition);
        }

        public Stream<MultiDeBruijnVertex> allVertices = kmerMap.values().stream()
                .flatMap(clusterer -> clusterer.getVertices().stream());

    }

    private static class PositionClusterer {
        private final Kmer kmer;

        // this will always get sorted from greatest to least count whenever thing get too wrong
        private final List<Pair<Integer, MutableInt>> clustersAndCounts = new ArrayList<>(1);

        private final int tolerance;

        private boolean finalized;

        private final List<MultiDeBruijnVertex> vertices = new ArrayList<>(1);


        public PositionClusterer(final int tolerance, final Kmer kmer) {
            this.tolerance = tolerance;
            this.kmer = kmer;
        }

        public List<MultiDeBruijnVertex> getVertices() { return vertices; }

        public void finish() {
            if (finalized) {
                return;
            }
            finalized = true;
            sort();
            clustersAndCounts.forEach(pair -> vertices.add(new MultiDeBruijnVertex(kmer.bases())));
        }

        public Optional<MultiDeBruijnVertex> getVertex(final int position) {
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

        public Optional<MultiDeBruijnVertex> getVertexFromMaxPosition(final int maxPosition) {
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

        public Optional<MultiDeBruijnVertex> getVertexFromMinPosition(final int minPosition) {
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



}
