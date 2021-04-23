package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;

import java.io.File;
import java.util.Collections;

public class MutectReadThreadingAssemblerArgumentCollection extends ReadThreadingAssemblerArgumentCollection {
    private static final long serialVersionUID = 5304L;

    /**
     * A single edge multiplicity cutoff for pruning doesn't work in samples with variable depths, for example exomes
     * and RNA.  This parameter disables the probabilistic algorithm for pruning the assembly graph that considers the
     * likelihood that each chain in the graph comes from real variation, and instead uses a simple multiplicity cutoff.
     */
    @Advanced
    @Argument(fullName="disable-adaptive-pruning", doc = "Disable the adaptive algorithm for pruning paths in the graph", optional = true)
    public boolean disableAdaptivePruning = false;

    @Override
    public ReadThreadingAssembler makeReadThreadingAssembler() {
        final ReadThreadingAssembler assemblyEngine = new ReadThreadingAssembler(maxNumHaplotypesInPopulation, Collections.unmodifiableList(kmerSizes),
                dontIncreaseKmerSizesForCycles, allowNonUniqueKmersInRef, numPruningSamples, disableAdaptivePruning ? minPruneFactor : 0,
                !disableAdaptivePruning, initialErrorRateForPruning, pruningLogOddsThreshold, pruningSeedingLogOddsThreshold, maxUnprunedVariants, useLinkedDeBruijnGraph,
                enableLegacyGraphCycleDetection, minMatchingBasesToDanglingEndRecovery);
        assemblyEngine.setDebugGraphTransformations(debugGraphTransformations);
        assemblyEngine.setRecoverDanglingBranches(true);
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
