package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;

import java.io.File;
import java.util.Collections;

public class HaplotypeCallerReadThreadingAssemblerArgumentCollection extends ReadThreadingAssemblerArgumentCollection {
    private static final long serialVersionUID = 6520834L;
    public static final String ADAPTIVE_PRUNING_LONG_NAME = "adaptive-pruning";
    public static final String DO_NOT_RECOVER_DANGLING_BRANCHES_LONG_NAME = "do-not-recover-dangling-branches";
    public static final String RECOVER_DANGLING_HEADS_LONG_NAME = "recover-dangling-heads";
    /**
     * A single edge multiplicity cutoff for pruning doesn't work in samples with variable depths, for example exomes
     * and RNA.  This parameter enables the probabilistic algorithm for pruning the assembly graph that considers the
     * likelihood that each chain in the graph comes from real variation.
     */
    @Advanced
    @Argument(fullName= ADAPTIVE_PRUNING_LONG_NAME, doc = "Use Mutect2's adaptive graph pruning algorithm", optional = true)
    public boolean useAdaptivePruning = false;

    /**
     * By default, the read threading assembler will attempt to recover dangling heads and tails. See the `minDanglingBranchLength` argument documentation for more details.
     */
    @Hidden
    @Argument(fullName= DO_NOT_RECOVER_DANGLING_BRANCHES_LONG_NAME, doc="Disable dangling head and tail recovery", optional = true)
    public boolean doNotRecoverDanglingBranches = false;

    /**
     * As of version 3.3, this argument is no longer needed because dangling end recovery is now the default behavior. See GATK 3.3 release notes for more details.
     */
    @Deprecated
    @Argument(fullName= RECOVER_DANGLING_HEADS_LONG_NAME, doc="This argument is deprecated since version 3.3", optional = true)
    public boolean DEPRECATED_RecoverDanglingHeads = false;


    @Override
    public ReadThreadingAssembler makeReadThreadingAssembler() {
        final ReadThreadingAssembler assemblyEngine = new ReadThreadingAssembler(maxNumHaplotypesInPopulation, Collections.unmodifiableList(kmerSizes),
                dontIncreaseKmerSizesForCycles, allowNonUniqueKmersInRef, numPruningSamples, useAdaptivePruning ? 0 : minPruneFactor,
                useAdaptivePruning, initialErrorRateForPruning, pruningLogOddsThreshold, pruningSeedingLogOddsThreshold, maxUnprunedVariants, useLinkedDeBruijnGraph,
                enableLegacyGraphCycleDetection, minMatchingBasesToDanglingEndRecovery);
        assemblyEngine.setDebugGraphTransformations(debugGraphTransformations);
        assemblyEngine.setRecoverDanglingBranches(!doNotRecoverDanglingBranches);
        assemblyEngine.setRecoverAllDanglingBranches(recoverAllDanglingBranches);
        assemblyEngine.setMinDanglingBranchLength(minDanglingBranchLength);
        assemblyEngine.setArtificialHaplotypeRecoveryMode(disableArtificialHaplotypeRecovery);

        if ( graphOutput != null ) {
            assemblyEngine.setGraphWriter(new File(graphOutput));
        }
        if ( haplotypeHistogramOutput != null ) {
            assemblyEngine.setDebugHistogramOutput(new File(haplotypeHistogramOutput));
        }

        return assemblyEngine;
    }
}
