package org.broadinstitute.hellbender.tools.genomicsdb;

import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;

import java.nio.file.Path;

/**
 * Encapsulates the GenomicsDB-specific options relevant to the FeatureDataSource
 */
public final class GenomicsDBOptions {
    final private Path reference;
    final private boolean callGenotypes;
    final private int maxDiploidAltAllelesThatCanBeGenotyped;
    final private int maxGenotypeCount;
    final private boolean useBCFCodec;
    final private boolean sharedPosixFSOptimizations;
    final private boolean useGcsHdfsConnector;

    public GenomicsDBOptions() {
        this(null);
    }

    public GenomicsDBOptions(final Path reference) {
        this(reference, new GenomicsDBArgumentCollection());
    }

    public GenomicsDBOptions(final Path reference, GenomicsDBArgumentCollection genomicsdbArgs) {
        this(reference, genomicsdbArgs, new GenotypeCalculationArgumentCollection(), false);
    }

    public GenomicsDBOptions(final Path reference, final GenomicsDBArgumentCollection genomicsdbArgs,
                             final GenotypeCalculationArgumentCollection genotypeCalcArgs, final boolean forceCallGenotypes) {
        this.reference = reference;
        this.callGenotypes = genomicsdbArgs.callGenotypes || forceCallGenotypes;
        this.maxDiploidAltAllelesThatCanBeGenotyped = genomicsdbArgs.maxDiploidAltAllelesThatCanBeGenotyped;
        this.maxGenotypeCount = genotypeCalcArgs.MAX_GENOTYPE_COUNT;
        this.useBCFCodec = genomicsdbArgs.useBCFCodec;
        this.sharedPosixFSOptimizations = genomicsdbArgs.sharedPosixFSOptimizations;
        this.useGcsHdfsConnector = genomicsdbArgs.useGcsHdfsConnector;
    }

    public Path getReference() {
        return reference;
    }

    public boolean doCallGenotypes() {
        return callGenotypes;
    }

    public int getMaxDiploidAltAllelesThatCanBeGenotyped() {
        return maxDiploidAltAllelesThatCanBeGenotyped;
    }

    public int getMaxGenotypeCount() {
        return maxGenotypeCount;
    }

    public boolean useBCFCodec() {
        return useBCFCodec;
    }

    public boolean sharedPosixFSOptimizations() {
        return sharedPosixFSOptimizations;
    }

    public boolean useGcsHdfsConnector() {
        return useGcsHdfsConnector;
    }
}
