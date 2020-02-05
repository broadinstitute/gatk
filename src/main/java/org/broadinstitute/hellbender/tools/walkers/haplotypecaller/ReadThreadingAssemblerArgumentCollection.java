package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.collect.Lists;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.io.Serializable;
import java.util.List;

/**
 * Set of arguments related to the {@link org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler}
 */
public abstract class ReadThreadingAssemblerArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final double DEFAULT_PRUNING_LOG_ODDS_THRESHOLD = MathUtils.log10ToLog(1.0);

    public static final String ERROR_CORRECT_READS_LONG_NAME = "error-correct-reads";

    public static final String CAPTURE_ASSEMBLY_FAILURE_BAM_LONG_NAME = "capture-assembly-failure-bam";
    public static final String KMER_SIZE_LONG_NAME = "kmer-size";
    public static final String DONT_INCREASE_KMER_SIZE_LONG_NAME = "dont-increase-kmer-sizes-for-cycles";

    // -----------------------------------------------------------------------------------------------
    // arguments to control internal behavior of the read threading assembler
    // -----------------------------------------------------------------------------------------------

    /**
     * Multiple kmer sizes can be specified, using e.g. `--kmer-size 10 --kmer-size 25`.
     */
    @Advanced
    @Argument(fullName= KMER_SIZE_LONG_NAME, doc="Kmer size to use in the read threading assembler", optional = true)
    public List<Integer> kmerSizes = Lists.newArrayList(10,25);

    /**
     * When graph cycles are detected, the normal behavior is to increase kmer sizes iteratively until the cycles are
     * resolved. Disabling this behavior may cause the program to give up on assembling the ActiveRegion.
     */
    @Advanced
    @Argument(fullName= DONT_INCREASE_KMER_SIZE_LONG_NAME, doc="Disable iterating over kmer sizes when graph cycles are detected", optional = true)
    public boolean dontIncreaseKmerSizesForCycles = false;

    /**
     * By default, the program does not allow processing of reference sections that contain non-unique kmers. Disabling
     * this check may cause problems in the assembly graph.
     */
    @Advanced
    @Argument(fullName="allow-non-unique-kmers-in-ref", doc="Allow graphs that have non-unique kmers in the reference", optional = true)
    public boolean allowNonUniqueKmersInRef = false;

    /**
     * If fewer samples than the specified number pass the minPruning threshold for a given path, that path will be eliminated from the graph.
     */
    @Advanced
    @Argument(fullName="num-pruning-samples", doc="Number of samples that must pass the minPruning threshold", optional = true)
    public int numPruningSamples = 1;

    /**
     * When constructing the assembly graph we are often left with "dangling" branches.  The assembly engine attempts to rescue these branches
     * by merging them back into the main graph.  This argument describes the minimum length of a dangling branch needed for the engine to
     * try to rescue it.  A smaller number here will lead to higher sensitivity to real variation but also to a higher number of false positives.
     */
    @Advanced
    @Argument(fullName="min-dangling-branch-length", doc="Minimum length of a dangling branch to attempt recovery", optional = true)
    public int minDanglingBranchLength = 4;

    /**
     * By default, the read threading assembler does not recover dangling branches that fork after splitting from the reference.  This argument
     * tells the assembly engine to recover all dangling branches.
     */
    @Advanced
    @Argument(fullName="recover-all-dangling-branches", doc="Recover all dangling branches", optional = true)
    public boolean recoverAllDanglingBranches = false;

    /**
     * The assembly graph can be quite complex, and could imply a very large number of possible haplotypes.  Each haplotype
     * considered requires N PairHMM evaluations if there are N reads across all samples.  In order to control the
     * run of the haplotype caller we only take maxNumHaplotypesInPopulation paths from the graph, in order of their
     * weights, no matter how many paths are possible to generate from the graph.  Putting this number too low
     * will result in dropping true variation because paths that include the real variant are not even considered.
     * You can consider increasing this number when calling organisms with high heterozygosity.
     */
    @Advanced
    @Argument(fullName="max-num-haplotypes-in-population", doc="Maximum number of haplotypes to consider for your population", optional = true)
    public int maxNumHaplotypesInPopulation = 128;

    /**
     * Paths with fewer supporting kmers than the specified threshold will be pruned from the graph.
     *
     * Be aware that this argument can dramatically affect the results of variant calling and should only be used with great caution.
     * Using a prune factor of 1 (or below) will prevent any pruning from the graph, which is generally not ideal; it can make the
     * calling much slower and even less accurate (because it can prevent effective merging of "tails" in the graph).  Higher values
     * tend to make the calling much faster, but also lowers the sensitivity of the results (because it ultimately requires higher
     * depth to produce calls).
     */
    @Advanced
    @Argument(fullName="min-pruning", doc = "Minimum support to not prune paths in the graph", optional = true)
    public int minPruneFactor = 2;

    /**
     * Initial base error rate guess for the probabilistic adaptive pruning model.  Results are not very sensitive to this
     * parameter because it is only a starting point from which the algorithm discovers the true error rate.
     */
    @Advanced
    @Argument(fullName="adaptive-pruning-initial-error-rate", doc = "Initial base error rate estimate for adaptive pruning", optional = true)
    public double initialErrorRateForPruning = 0.001;

    /**
     * Log-10 likelihood ratio threshold for adaptive pruning algorithm.
     */
    @Advanced
    @Argument(fullName="pruning-lod-threshold", doc = "Ln likelihood ratio threshold for adaptive pruning algorithm", optional = true)
    public double pruningLogOddsThreshold = DEFAULT_PRUNING_LOG_ODDS_THRESHOLD;

    /**
     * The maximum number of variants in graph the adaptive pruner will allow
     */
    @Advanced
    @Argument(fullName="max-unpruned-variants", doc = "Maximum number of variants in graph the adaptive pruner will allow", optional = true)
    public int maxUnprunedVariants = 100;

    /**
     * Disables graph simplification into a seq graph, opts to construct a proper De Brujin graph with potential loops
     *
     * NOTE: --linked-de-bruijn-graph is currently an experimental feature that does not directly match with
     *        the regular HaplotypeCaller. Specifically the haplotype finding code does not perform correctly at complicated
     *        sites. Use this mode at your own risk.
     */
    @Hidden
    @Argument(fullName="linked-de-bruijn-graph", doc = "If enabled, the Assembly Engine will construct a Linked De Brujin graph to recover better haplotypes", optional = true)
    public boolean useLinkedDeBrujinGraph = false;

    @Advanced
    @Argument(fullName="debug-assembly", shortName="debug", doc="Print out verbose debug information about each assembly region", optional = true)
    public boolean debugAssembly;

    @Hidden
    @Argument(fullName="debug-graph-transformations", doc="Write DOT formatted graph files out of the assembler for only this graph size", optional = true)
    public boolean debugGraphTransformations = false;

    /**
     * This argument is meant for debugging and is not immediately useful for normal analysis use.
     */
    @Argument(fullName="graph-output", shortName="graph", doc="Write debug assembly graph information to this file", optional = true)
    public String graphOutput = null;

    /**
     * This argument is meant for debugging and is not immediately useful for normal analysis use.
     */
    @Hidden
    @Argument(fullName="haplotype-debug-histogram-output", doc="Write debug assembly graph information to this file", optional = true)
    public String haplotypeHistogramOutput = null;

    @Hidden
    @Argument(fullName = CAPTURE_ASSEMBLY_FAILURE_BAM_LONG_NAME, doc = "Write a BAM called assemblyFailure.bam capturing all of the reads that were in the active region when the assembler failed for any reason", optional = true)
    public boolean captureAssemblyFailureBAM = false;

    //---------------------------------------------------------------------------------------------------------------
    //
    // Read Error Corrector Related Parameters
    //
    // ---------------------------------------------------------------------------------------------------------------
    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    @Hidden
    @Argument(fullName = ERROR_CORRECT_READS_LONG_NAME, doc = "Use an exploratory algorithm to error correct the kmers used during assembly", optional = true)
    public boolean errorCorrectReads = false;

    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    @Hidden
    @Argument(fullName="kmer-length-for-read-error-correction", doc = "Use an exploratory algorithm to error correct the kmers used during assembly", optional = true)
    public int kmerLengthForReadErrorCorrection = 25;

    @Hidden
    @Argument(fullName="min-observations-for-kmer-to-be-solid", doc = "A k-mer must be seen at least these times for it considered to be solid", optional = true)
    public int minObservationsForKmerToBeSolid = 20;

    public abstract ReadThreadingAssembler makeReadThreadingAssembler();
}
