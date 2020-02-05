package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResult;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadErrorCorrector;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.*;
import org.broadinstitute.hellbender.utils.Histogram;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public final class ReadThreadingAssembler {
    private static final Logger logger = LogManager.getLogger(ReadThreadingAssembler.class);

    private static final int DEFAULT_NUM_PATHS_PER_GRAPH = 128;
    private static final int KMER_SIZE_ITERATION_INCREASE = 10;
    private static final int MAX_KMER_ITERATIONS_TO_ATTEMPT = 6;

    /** The min and max kmer sizes to try when building the graph. */
    private final List<Integer> kmerSizes;
    private final boolean dontIncreaseKmerSizesForCycles;
    private final boolean allowNonUniqueKmersInRef;
    private final boolean generateSeqGraph;
    private final int numPruningSamples;
    private final int numBestHaplotypesPerGraph;

    private boolean removePathsNotConnectedToRef = true;
    private boolean justReturnRawGraph = false;

    /**
     * If false, we will only write out a region around the reference source
     */
    private static final boolean PRINT_FULL_GRAPH_FOR_DEBUGGING = true;
    private static final byte DEFAULT_MIN_BASE_QUALITY_TO_USE = (byte) 10;
    private static final int MIN_HAPLOTYPE_REFERENCE_LENGTH = 30;

    private boolean debug = false;
    private boolean debugGraphTransformations = false;
    private boolean recoverDanglingBranches = true;
    private boolean recoverAllDanglingBranches = false;
    private int minDanglingBranchLength = 0;
    
    protected byte minBaseQualityToUseInAssembly = DEFAULT_MIN_BASE_QUALITY_TO_USE;
    private int pruneFactor;
    private final ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> chainPruner;

    private File debugGraphOutputPath = null;  //Where to write debug graphs, if unset it defaults to the current working dir
    private File graphOutputPath = null;
    private File graphHaplotypeHistogramPath = null;
    private Histogram haplotypeHistogram = null;
    private Histogram kmersUsedHistogram = null;

    public ReadThreadingAssembler(final int maxAllowedPathsForReadThreadingAssembler, final List<Integer> kmerSizes,
                                  final boolean dontIncreaseKmerSizesForCycles, final boolean allowNonUniqueKmersInRef,
                                  final int numPruningSamples, final int pruneFactor, final boolean useAdaptivePruning,
                                  final double initialErrorRateForPruning, final double pruningLogOddsThreshold,
                                  final int maxUnprunedVariants, final boolean useLinkedDebrujinGraphs) {
        Utils.validateArg( maxAllowedPathsForReadThreadingAssembler >= 1, "numBestHaplotypesPerGraph should be >= 1 but got " + maxAllowedPathsForReadThreadingAssembler);
        this.kmerSizes = kmerSizes;
        this.dontIncreaseKmerSizesForCycles = dontIncreaseKmerSizesForCycles;
        this.allowNonUniqueKmersInRef = allowNonUniqueKmersInRef;
        this.numPruningSamples = numPruningSamples;
        this.pruneFactor = pruneFactor;
        this.generateSeqGraph = !useLinkedDebrujinGraphs;
        chainPruner = useAdaptivePruning ? new AdaptiveChainPruner<>(initialErrorRateForPruning, pruningLogOddsThreshold, maxUnprunedVariants) :
                new LowWeightChainPruner<>(pruneFactor);
        numBestHaplotypesPerGraph = maxAllowedPathsForReadThreadingAssembler;
    }

    @VisibleForTesting
    ReadThreadingAssembler(final int maxAllowedPathsForReadThreadingAssembler, final List<Integer> kmerSizes, final int pruneFactor) {
        this(maxAllowedPathsForReadThreadingAssembler, kmerSizes, true, true, 1, pruneFactor, false, 0.001, 2, Integer.MAX_VALUE, false);
    }

    @VisibleForTesting
    ReadThreadingAssembler() {
        this(DEFAULT_NUM_PATHS_PER_GRAPH, Arrays.asList(25), 2);
    }

    /**
     * Main entry point into the assembly engine. Build a set of deBruijn graphs out of the provided reference sequence and list of reads
     * @param assemblyRegion              AssemblyRegion object holding the reads which are to be used during assembly
     * @param refHaplotype              reference haplotype object
     * @param fullReferenceWithPadding  byte array holding the reference sequence with padding
     * @param refLoc                    GenomeLoc object corresponding to the reference sequence with padding
     * @param readErrorCorrector        a ReadErrorCorrector object, if read are to be corrected before assembly. Can be null if no error corrector is to be used.
     * @param aligner                   {@link SmithWatermanAligner} used to align dangling ends in assembly graphs to the reference sequence
     * @return                          the resulting assembly-result-set
     */
    public AssemblyResultSet runLocalAssembly(final AssemblyRegion assemblyRegion,
                                              final Haplotype refHaplotype,
                                              final byte[] fullReferenceWithPadding,
                                              final SimpleInterval refLoc,
                                              final ReadErrorCorrector readErrorCorrector,
                                              final SAMFileHeader header,
                                              final SmithWatermanAligner aligner) {
        Utils.nonNull(assemblyRegion, "Assembly engine cannot be used with a null AssemblyRegion.");
        Utils.nonNull(assemblyRegion.getPaddedSpan(), "Active region must have an extended location.");
        Utils.nonNull(refHaplotype, "Reference haplotype cannot be null.");
        Utils.nonNull(fullReferenceWithPadding, "fullReferenceWithPadding");
        Utils.nonNull(refLoc, "refLoc");
        Utils.nonNull(aligner, "aligner");
        Utils.validateArg( fullReferenceWithPadding.length == refLoc.size(), "Reference bases and reference loc must be the same size.");
        ParamUtils.isPositiveOrZero(pruneFactor, "Pruning factor cannot be negative");

        // error-correct reads before clipping low-quality tails: some low quality bases might be good and we want to recover them
        final List<GATKRead> correctedReads;
        if ( readErrorCorrector != null ) {
            // now correct all reads in active region after filtering/downsampling
            // Note that original reads in active region are NOT modified by default, since they will be used later for GL computation,
            // and we only want the read-error corrected reads for graph building.
            readErrorCorrector.addReadsToKmers(assemblyRegion.getReads());
            correctedReads = new ArrayList<>(readErrorCorrector.correctReads(assemblyRegion.getReads()));
        } else {
            correctedReads = assemblyRegion.getReads();
        }

        final List<AbstractReadThreadingGraph> nonRefRTGraphs = new LinkedList<>();
        final List<SeqGraph> nonRefSeqGraphs = new LinkedList<>();
        final AssemblyResultSet resultSet = new AssemblyResultSet();
        resultSet.setRegionForGenotyping(assemblyRegion);
        resultSet.setFullReferenceWithPadding(fullReferenceWithPadding);
        resultSet.setPaddedReferenceLoc(refLoc);
        final SimpleInterval activeRegionExtendedLocation = assemblyRegion.getPaddedSpan();
        refHaplotype.setGenomeLocation(activeRegionExtendedLocation);
        resultSet.add(refHaplotype);
        final Map<AbstractReadThreadingGraph,AssemblyResult> assemblyResultByRTGraph = new HashMap<>();
        final Map<SeqGraph,AssemblyResult> assemblyResultBySeqGraph = new HashMap<>();
        // create the graphs by calling our subclass assemble method
        for ( final AssemblyResult result : assemble(correctedReads, refHaplotype, header, aligner) ) {
            if ( result.getStatus() == AssemblyResult.Status.ASSEMBLED_SOME_VARIATION ) {
                // do some QC on the graph
                if (generateSeqGraph) {
                    sanityCheckGraph(result.getSeqGraph(), refHaplotype);
                    // add it to graphs with meaningful non-reference features
                    assemblyResultBySeqGraph.put(result.getSeqGraph(),result);
                    nonRefSeqGraphs.add(result.getSeqGraph());
                } else {
                    sanityCheckGraph(result.getThreadingGraph(), refHaplotype);
                    result.getThreadingGraph().postProcessForHaplotypeFinding(debugGraphOutputPath, refHaplotype.getLocation());
                    // add it to graphs with meaningful non-reference features
                    assemblyResultByRTGraph.put(result.getThreadingGraph(),result);
                    nonRefRTGraphs.add(result.getThreadingGraph());
                }

                if (graphHaplotypeHistogramPath != null) {
                    kmersUsedHistogram.add((double)result.getKmerSize());
                }
            }
        }

        // add assembled alt haplotypes to the {@code resultSet}
        if (generateSeqGraph) {
            findBestPaths(nonRefSeqGraphs, refHaplotype, refLoc, activeRegionExtendedLocation, assemblyResultBySeqGraph, resultSet, aligner);
        } else {
            findBestPaths(nonRefRTGraphs, refHaplotype, refLoc, activeRegionExtendedLocation, assemblyResultByRTGraph, resultSet, aligner);
        }

        // print the graphs if the appropriate debug option has been turned on
        if ( graphOutputPath != null ) { printGraphs(nonRefSeqGraphs); }
        if ( graphHaplotypeHistogramPath != null ) { haplotypeHistogram.add((double)resultSet.getHaplotypeCount()); }

        return resultSet;
    }

    private <V extends  BaseVertex, E extends BaseEdge, T extends BaseGraph<V, E>>
    List<Haplotype> findBestPaths(final Collection<T> graphs, final Haplotype refHaplotype, final SimpleInterval refLoc, final SimpleInterval activeRegionWindow,
                                          final Map<T, AssemblyResult> assemblyResultByGraph, final AssemblyResultSet assemblyResultSet, final SmithWatermanAligner aligner) {
        // add the reference haplotype separately from all the others to ensure that it is present in the list of haplotypes
        final Set<Haplotype> returnHaplotypes = new LinkedHashSet<>();

        final int activeRegionStart = refHaplotype.getAlignmentStartHapwrtRef();
        int failedCigars = 0;

        for( final BaseGraph<V, E> graph : graphs ) {
            final V source = graph.getReferenceSourceVertex();
            final V sink = graph.getReferenceSinkVertex();
            Utils.validateArg( source != null && sink != null, () -> "Both source and sink cannot be null but got " + source + " and sink " + sink + " for graph " + graph);

            for (final KBestHaplotype<V, E> kBestHaplotype :
                    (generateSeqGraph ?
                            new GraphBasedKBestHaplotypeFinder<>(graph,source,sink) :
                            new JunctionTreeKBestHaplotypeFinder<>(graph,source,sink, JunctionTreeKBestHaplotypeFinder.DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE))
                            .findBestHaplotypes(numBestHaplotypesPerGraph)) {
                final Haplotype h = kBestHaplotype.haplotype();
                if( !returnHaplotypes.contains(h) ) {
                    // TODO this score seems to be irrelevant at this point...
                    if (kBestHaplotype.isReference()) {
                        refHaplotype.setScore(kBestHaplotype.score());
                    }
                    final Cigar cigar = CigarUtils.calculateCigar(refHaplotype.getBases(), h.getBases(), aligner, SWOverhangStrategy.SOFTCLIP);

                    if ( cigar == null ) {
                        failedCigars++; // couldn't produce a meaningful alignment of haplotype to reference, fail quietly
                        continue;
                    } else if( cigar.isEmpty() ) {
                        throw new IllegalStateException("Smith-Waterman alignment failure. Cigar = " + cigar + " with reference length " + cigar.getReferenceLength() +
                                " but expecting reference length of " + refHaplotype.getCigar().getReferenceLength());
                    } else if ( pathIsTooDivergentFromReference(cigar) || cigar.getReferenceLength() < MIN_HAPLOTYPE_REFERENCE_LENGTH ) {
                        // N cigar elements means that a bubble was too divergent from the reference so skip over this path
                        continue;
                    } else if( cigar.getReferenceLength() != refHaplotype.getCigar().getReferenceLength() ) { // SW failure
                        throw new IllegalStateException("Smith-Waterman alignment failure. Cigar = " + cigar + " with reference length "
                                + cigar.getReferenceLength() + " but expecting reference length of " + refHaplotype.getCigar().getReferenceLength()
                                + " ref = " + refHaplotype + " path " + new String(h.getBases()));
                    }

                    h.setCigar(cigar);
                    h.setAlignmentStartHapwrtRef(activeRegionStart);
                    h.setGenomeLocation(activeRegionWindow);
                    returnHaplotypes.add(h);
                    assemblyResultSet.add(h, assemblyResultByGraph.get(graph));

                    if ( debug ) {
                        logger.info("Adding haplotype " + h.getCigar() + " from graph with kmer " + graph.getKmerSize());
                    }
                }
            }
        }

        // Make sure that the ref haplotype is amongst the return haplotypes and calculate its score as
        // the first returned by any finder.
        if (!returnHaplotypes.contains(refHaplotype)) {
            returnHaplotypes.add(refHaplotype);
        }

        if (failedCigars != 0) {
            logger.debug(String.format("failed to align some haplotypes (%d) back to the reference (loc=%s); these will be ignored.", failedCigars, refLoc.toString()));
        }

        if ( debug ) {
            if( returnHaplotypes.size() > 1 ) {
                logger.info("Found " + returnHaplotypes.size() + " candidate haplotypes of " + returnHaplotypes.size() + " possible combinations to evaluate every read against.");
            } else {
                logger.info("Found only the reference haplotype in the assembly graph.");
            }
            for( final Haplotype h : returnHaplotypes ) {
                logger.info( h.toString() );
                logger.info( "> Cigar = " + h.getCigar() + " : " + h.getCigar().getReferenceLength() + " score " + h.getScore() + " ref " + h.isReference());
            }
        }

        return new ArrayList<>(returnHaplotypes);

    }
    /**
     * We use CigarOperator.N as the signal that an incomplete or too divergent bubble was found during bubble traversal
     * @param c the cigar to test
     * @return  true if we should skip over this path
     */
    private static boolean pathIsTooDivergentFromReference(final Cigar c) {
        return c.getCigarElements().stream().anyMatch(ce -> ce.getOperator() == CigarOperator.N);
    }

    /**
     * Print graph to file NOTE this requires that debugGraphTransformations be enabled.
     *
     * @param graph the graph to print
     * @param fileName the name to give the graph file
     */
    private void printDebugGraphTransform(final BaseGraph<?,?> graph, final String fileName) {
        File outputFile = new File(debugGraphOutputPath, fileName);
        if ( PRINT_FULL_GRAPH_FOR_DEBUGGING ) {
            graph.printGraph(outputFile, pruneFactor);
        } else {
            graph.subsetToRefSource(10).printGraph(outputFile, pruneFactor);
        }
    }

    private AssemblyResult getResultSetForRTGraph(final AbstractReadThreadingGraph rtGraph) {

        // The graph has degenerated in some way, so the reference source and/or sink cannot be id'd.  Can
        // happen in cases where for example the reference somehow manages to acquire a cycle, or
        // where the entire assembly collapses back into the reference sequence.
        if ( rtGraph.getReferenceSourceVertex() == null || rtGraph.getReferenceSinkVertex() == null ) {
            return new AssemblyResult(AssemblyResult.Status.JUST_ASSEMBLED_REFERENCE, null, rtGraph);
        }

        rtGraph.removePathsNotConnectedToRef();

        return new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION, null, rtGraph);
    }

    // Performs the various transformations necessary on a sequence graph
    private AssemblyResult cleanupSeqGraph(final SeqGraph seqGraph, final Haplotype refHaplotype) {
        if (debugGraphTransformations) {
            printDebugGraphTransform(seqGraph, refHaplotype.getLocation() + "-sequenceGraph."+seqGraph.getKmerSize()+".1.0.non_ref_removed.dot");
        }

        // the very first thing we need to do is zip up the graph, or pruneGraph will be too aggressive
        seqGraph.zipLinearChains();
        if (debugGraphTransformations) {
            printDebugGraphTransform(seqGraph, refHaplotype.getLocation() + "-sequenceGraph."+seqGraph.getKmerSize()+".1.1.zipped.dot");
        }

        // now go through and prune the graph, removing vertices no longer connected to the reference chain
        seqGraph.removeSingletonOrphanVertices();
        seqGraph.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();

        if (debugGraphTransformations) {
            printDebugGraphTransform(seqGraph, refHaplotype.getLocation() + "-sequenceGraph."+seqGraph.getKmerSize()+".1.2.pruned.dot");
        }
        seqGraph.simplifyGraph();
        if (debugGraphTransformations) {
            printDebugGraphTransform(seqGraph, refHaplotype.getLocation() + "-sequenceGraph."+seqGraph.getKmerSize()+".1.3.merged.dot");
        }

        // The graph has degenerated in some way, so the reference source and/or sink cannot be id'd.  Can
        // happen in cases where for example the reference somehow manages to acquire a cycle, or
        // where the entire assembly collapses back into the reference sequence.
        if ( seqGraph.getReferenceSourceVertex() == null || seqGraph.getReferenceSinkVertex() == null ) {
            return new AssemblyResult(AssemblyResult.Status.JUST_ASSEMBLED_REFERENCE, seqGraph, null);
        }

        seqGraph.removePathsNotConnectedToRef();
        seqGraph.simplifyGraph();
        if ( seqGraph.vertexSet().size() == 1 ) {
            // we've perfectly assembled into a single reference haplotype, add a empty seq vertex to stop
            // the code from blowing up.
            // TODO -- ref properties should really be on the vertices, not the graph itself
            final SeqVertex complete = seqGraph.vertexSet().iterator().next();
            final SeqVertex dummy = new SeqVertex("");
            seqGraph.addVertex(dummy);
            seqGraph.addEdge(complete, dummy, new BaseEdge(true, 0));
        }
        if (debugGraphTransformations) {
            printDebugGraphTransform(seqGraph, refHaplotype.getLocation() + "-sequenceGraph."+seqGraph.getKmerSize()+".1.4.final.dot");
        }
        return new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION, seqGraph, null);
    }

    /**
     * Perform general QC on the graph to make sure something hasn't gone wrong during assembly
     * @param graph the graph to check
     * @param refHaplotype the reference haplotype
     */
    private static <T extends BaseVertex, E extends BaseEdge> void sanityCheckGraph(final BaseGraph<T,E> graph, final Haplotype refHaplotype) {
        sanityCheckReferenceGraph(graph, refHaplotype);
    }

    /**
     * Make sure the reference sequence is properly represented in the provided graph
     *
     * @param graph the graph to check
     * @param refHaplotype the reference haplotype
     */
    private static <T extends BaseVertex, E extends BaseEdge> void sanityCheckReferenceGraph(final BaseGraph<T,E> graph, final Haplotype refHaplotype) {
        if( graph.getReferenceSourceVertex() == null ) {
            throw new IllegalStateException("All reference graphs must have a reference source vertex.");
        }
        if( graph.getReferenceSinkVertex() == null ) {
            throw new IllegalStateException("All reference graphs must have a reference sink vertex.");
        }
        if( !Arrays.equals(graph.getReferenceBytes(graph.getReferenceSourceVertex(), graph.getReferenceSinkVertex(), true, true), refHaplotype.getBases()) ) {
            throw new IllegalStateException("Mismatch between the reference haplotype and the reference assembly graph path. for graph " + graph +
                    " graph = " + new String(graph.getReferenceBytes(graph.getReferenceSourceVertex(), graph.getReferenceSinkVertex(), true, true)) +
                    " haplotype = " + new String(refHaplotype.getBases())
            );
        }
    }

    private static void addResult(final Collection<AssemblyResult> results, final AssemblyResult maybeNullResult) {
        if ( maybeNullResult != null ) {
            results.add(maybeNullResult);
        }
    }

    /**
     * Given reads and a reference haplotype give us graphs to use for constructing
     * non-reference haplotypes.
     *
     * @param reads the reads we're going to assemble
     * @param refHaplotype the reference haplotype
     * @param aligner {@link SmithWatermanAligner} used to align dangling ends in assembly graphs to the reference sequence
     * @return a non-null list of reads
     */
    @VisibleForTesting
    List<AssemblyResult> assemble(final List<GATKRead> reads, final Haplotype refHaplotype, final SAMFileHeader header, final SmithWatermanAligner aligner) {
        final List<AssemblyResult> results = new LinkedList<>();

        // first, try using the requested kmer sizes
        for ( final int kmerSize : kmerSizes ) {
            addResult(results, createGraph(reads, refHaplotype, kmerSize, dontIncreaseKmerSizesForCycles, allowNonUniqueKmersInRef, header, aligner));
        }

        // if none of those worked, iterate over larger sizes if allowed to do so
        if ( results.isEmpty() && !dontIncreaseKmerSizesForCycles ) {
            int kmerSize = arrayMaxInt(kmerSizes) + KMER_SIZE_ITERATION_INCREASE;
            int numIterations = 1;
            while ( results.isEmpty() && numIterations <= MAX_KMER_ITERATIONS_TO_ATTEMPT ) {
                // on the last attempt we will allow low complexity graphs
                final boolean lastAttempt = numIterations == MAX_KMER_ITERATIONS_TO_ATTEMPT;
                addResult(results, createGraph(reads, refHaplotype, kmerSize, lastAttempt, lastAttempt, header, aligner));
                kmerSize += KMER_SIZE_ITERATION_INCREASE;
                numIterations++;
            }
        }

        return results;
    }

    private static int arrayMaxInt(final List<Integer> array) {
        return array.stream().mapToInt(Integer::intValue).max().orElseThrow(() -> new IllegalArgumentException("Array size cannot be 0!"));
    }

    /**
     * Creates the sequence graph for the given kmerSize
     *
     * @param reads            reads to use
     * @param refHaplotype     reference haplotype
     * @param kmerSize         kmer size
     * @param allowLowComplexityGraphs if true, do not check for low-complexity graphs
     * @param allowNonUniqueKmersInRef if true, do not fail if the reference has non-unique kmers
     * @param aligner {@link SmithWatermanAligner} used to align dangling ends to the reference sequence
     * @return sequence graph or null if one could not be created (e.g. because it contains cycles or too many paths or is low complexity)
     */
    private AssemblyResult createGraph(final Iterable<GATKRead> reads,
                                       final Haplotype refHaplotype,
                                       final int kmerSize,
                                       final boolean allowLowComplexityGraphs,
                                       final boolean allowNonUniqueKmersInRef,
                                       final SAMFileHeader header,
                                       final SmithWatermanAligner aligner) {
        if ( refHaplotype.length() < kmerSize ) {
            // happens in cases where the assembled region is just too small
            return new AssemblyResult(AssemblyResult.Status.FAILED, null, null);
        }

        if ( !allowNonUniqueKmersInRef && !ReadThreadingGraph.determineNonUniqueKmers(new ReadThreadingGraph.SequenceForKmers("ref", refHaplotype.getBases(), 0, refHaplotype.getBases().length, 1, true), kmerSize).isEmpty() ) {
            if ( debug ) {
                logger.info("Not using kmer size of " + kmerSize + " in read threading assembler because reference contains non-unique kmers");
            }
            return null;
        }

        // TODO figure out how you want to hook this in
        final AbstractReadThreadingGraph rtgraph = generateSeqGraph ? new ReadThreadingGraph(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly, numPruningSamples) :
                new JunctionTreeLinkedDeBruinGraph(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly, numPruningSamples);

        rtgraph.setThreadingStartOnlyAtExistingVertex(!recoverDanglingBranches);

        // add the reference sequence to the graph
        rtgraph.addSequence("ref", refHaplotype.getBases(), true);

        // Next pull kmers out of every read and throw them on the graph
        for( final GATKRead read : reads ) {
            rtgraph.addRead(read, header);
        }

        // actually build the read threading graph
        rtgraph.buildGraphIfNecessary();

        // It's important to prune before recovering dangling ends so that we don't waste time recovering bad ends.
        // It's also important to prune before checking for cycles so that sequencing errors don't create false cycles
        // and unnecessarily abort assembly
        chainPruner.pruneLowWeightChains(rtgraph);

        // sanity check: make sure there are no cycles in the graph, unless we are in experimental mode
        if ( generateSeqGraph && rtgraph.hasCycles() ) {
            if ( debug ) {
                logger.info("Not using kmer size of " + kmerSize + " in read threading assembler because it contains a cycle");
            }
            return null;
        }

        // sanity check: make sure the graph had enough complexity with the given kmer
        if ( ! allowLowComplexityGraphs && rtgraph.isLowQualityGraph() ) {
            if ( debug ) {
                logger.info("Not using kmer size of " + kmerSize + " in read threading assembler because it does not produce a graph with enough complexity");
            }
            return null;
        }

        final AssemblyResult result = getAssemblyResult(refHaplotype, kmerSize, rtgraph, aligner);
        // check whether recovering dangling ends created cycles
        if (recoverAllDanglingBranches && rtgraph.hasCycles()) {
            return null;
        }
        return result;
    }

    private AssemblyResult getAssemblyResult(final Haplotype refHaplotype, final int kmerSize, final AbstractReadThreadingGraph rtgraph, final SmithWatermanAligner aligner) {
        if (debugGraphTransformations) {
            printDebugGraphTransform(rtgraph, refHaplotype.getLocation() + "-sequenceGraph." + kmerSize + ".0.0.raw_readthreading_graph.dot");
        }

        // look at all chains in the graph that terminate in a non-ref node (dangling sources and sinks) and see if
        // we can recover them by merging some N bases from the chain back into the reference
        if ( recoverDanglingBranches ) {
            rtgraph.recoverDanglingTails(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches, aligner);
            rtgraph.recoverDanglingHeads(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches, aligner);
        }

        // remove all heading and trailing paths
        if ( removePathsNotConnectedToRef ) {
            rtgraph.removePathsNotConnectedToRef();
        }

        if (debugGraphTransformations) {
            printDebugGraphTransform(rtgraph, refHaplotype.getLocation() + "-sequenceGraph." + kmerSize + ".0.1.cleaned_readthreading_graph.dot");
        }

        // Either return an assembly result with a sequence graph or with an unchanged sequence graph deptending on the kmer duplication behavior
        if (generateSeqGraph) {
            final SeqGraph initialSeqGraph = rtgraph.toSequenceGraph();
            if (debugGraphTransformations) {
                rtgraph.printGraph(new File(debugGraphOutputPath, refHaplotype.getLocation() + "-sequenceGraph." + kmerSize + ".0.1.initial_seqgraph.dot"), 10000);
            }

            // if the unit tests don't want us to cleanup the graph, just return the raw sequence graph
            if (justReturnRawGraph) {
                return new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION, initialSeqGraph, null);
            }

            if (debug) {
                logger.info("Using kmer size of " + rtgraph.getKmerSize() + " in read threading assembler");
            }
            if (debugGraphTransformations) {
                printDebugGraphTransform(initialSeqGraph, refHaplotype.getLocation() + "-sequenceGraph." + kmerSize + ".0.2.initial_seqgraph.dot");
            }
            initialSeqGraph.cleanNonRefPaths(); // TODO -- I don't this is possible by construction

            final AssemblyResult cleaned = cleanupSeqGraph(initialSeqGraph, refHaplotype);
            final AssemblyResult.Status status = cleaned.getStatus();
            return new AssemblyResult(status, cleaned.getSeqGraph(), rtgraph);

        } else {
            if (debugGraphTransformations) {
                rtgraph.printGraph(new File(debugGraphOutputPath, refHaplotype.getLocation() + "-sequenceGraph." + kmerSize + ".0.1.initial_seqgraph.dot"), 10000);
            }

            // if the unit tests don't want us to cleanup the graph, just return the raw sequence graph
            if (justReturnRawGraph) {
                return new AssemblyResult(AssemblyResult.Status.ASSEMBLED_SOME_VARIATION, null, rtgraph);
            }

            if (debug) {
                logger.info("Using kmer size of " + rtgraph.getKmerSize() + " in read threading assembler");
            }
            printDebugGraphTransform(rtgraph, refHaplotype.getLocation() + "-sequenceGraph." + kmerSize + ".0.2.initial_seqgraph.dot");
            rtgraph.cleanNonRefPaths(); // TODO -- I don't this is possible by construction

            final AssemblyResult cleaned = getResultSetForRTGraph(rtgraph);
            final AssemblyResult.Status status = cleaned.getStatus();
            return new AssemblyResult(status, null , cleaned.getThreadingGraph());
        }
    }

    @Override
    public String toString() {
        return "ReadThreadingAssembler{kmerSizes=" + kmerSizes + '}';
    }

    /**
     * Print the generated graphs to the graphWriter
     * @param graphs a non-null list of graphs to print out
     */
    private <T extends BaseGraph<?, ?>> void printGraphs(final List<T> graphs) {
        final int writeFirstGraphWithSizeSmallerThan = 50;

        try ( final PrintStream graphWriter = new PrintStream(graphOutputPath) ) {
            graphWriter.println("digraph assemblyGraphs {");
            for ( final T graph : graphs ) {
                if ( debugGraphTransformations && graph.getKmerSize() >= writeFirstGraphWithSizeSmallerThan ) {
                    logger.info("Skipping writing of graph with kmersize " + graph.getKmerSize());
                    continue;
                }

                graph.printGraph(graphWriter, false, pruneFactor);

                if ( debugGraphTransformations )
                    break;
            }

            graphWriter.println("}");
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile(graphOutputPath, e);
        }
    }

    /**
     * Print the generated graphs to the graphWriter
     */
    public void printDebugHistograms() {
        if (graphHaplotypeHistogramPath != null) {

            try (final PrintStream histogramWriter = new PrintStream(graphHaplotypeHistogramPath)) {
                histogramWriter.println("Histogram over the number of haplotypes recovered per active region:");
                histogramWriter.println(haplotypeHistogram.toString());

                histogramWriter.println("\nHistogram over the average kmer size used for assembly:");
                histogramWriter.println(kmersUsedHistogram.toString());

            } catch (IOException e) {
                throw new UserException.CouldNotCreateOutputFile(graphHaplotypeHistogramPath, e);
            }
        }
    }

    // -----------------------------------------------------------------------------------------------
    //
    // getter / setter routines for generic assembler properties
    //
    // -----------------------------------------------------------------------------------------------

    public void setGraphWriter(File graphOutputPath) {
        this.graphOutputPath = graphOutputPath;
    }

    public byte getMinBaseQualityToUseInAssembly() {
        return minBaseQualityToUseInAssembly;
    }

    public void setMinBaseQualityToUseInAssembly(byte minBaseQualityToUseInAssembly) {
        this.minBaseQualityToUseInAssembly = minBaseQualityToUseInAssembly;
    }

    public boolean isDebug() {
        return debug;
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    public boolean isDebugGraphTransformations() {
        return debugGraphTransformations;
    }

    public boolean isRecoverDanglingBranches() { return recoverDanglingBranches; }

    public void setDebugHistogramOutput(final File file) {
        this.graphHaplotypeHistogramPath = file;
        this.haplotypeHistogram = new Histogram(1.0);
        this.kmersUsedHistogram = new Histogram(1.0);
    }

    public void setDebugGraphTransformations(final boolean debugHaplotypeFinding) {
        this.debugGraphTransformations = debugHaplotypeFinding;
    }

    /**
     * Set where to write debug graph files if {@link ReadThreadingAssembler#debugGraphTransformations} == true
     */
    public void setDebugGraphOutputPath(File debugGraphOutputPath) {
        this.debugGraphOutputPath = debugGraphOutputPath;
    }

    public void setRecoverDanglingBranches(final boolean recoverDanglingBranches) {
        this.recoverDanglingBranches = recoverDanglingBranches;
    }

    public void setRecoverAllDanglingBranches(final boolean recoverAllDanglingBranches) {
        this.recoverAllDanglingBranches = recoverAllDanglingBranches;
        recoverDanglingBranches = true;
    }

    public void setMinDanglingBranchLength( final int minDanglingBranchLength ) {
        this.minDanglingBranchLength = minDanglingBranchLength;
    }

    @VisibleForTesting
    void setJustReturnRawGraph(final boolean justReturnRawGraph) {
        this.justReturnRawGraph = justReturnRawGraph;
    }

    public void setRemovePathsNotConnectedToRef(final boolean removePathsNotConnectedToRef) {
        this.removePathsNotConnectedToRef = removePathsNotConnectedToRef;
    }
}