package org.broadinstitute.hellbender.tools.genomicsdb;

import org.broadinstitute.hellbender.exceptions.UserException;
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
        this(reference, genomicsdbArgs, new GenotypeCalculationArgumentCollection());
    }

    public GenomicsDBOptions(final Path reference, final GenomicsDBArgumentCollection genomicsdbArgs,
                             final GenotypeCalculationArgumentCollection genotypeCalcArgs) {
        this.reference = reference;
        this.callGenotypes = genomicsdbArgs.callGenotypes;
        this.useBCFCodec = genomicsdbArgs.useBCFCodec;
        this.sharedPosixFSOptimizations = genomicsdbArgs.sharedPosixFSOptimizations;
        this.useGcsHdfsConnector = genomicsdbArgs.useGcsHdfsConnector;
        if (genomicsdbArgs.maxDiploidAltAllelesThatCanBeGenotyped - 1 < genotypeCalcArgs.maxAlternateAlleles) {  //-1 for <NON_REF>
            throw new UserException.BadInput("GenomicsDB max alternate alleles (" + GenomicsDBArgumentCollection.MAX_ALTS_LONG_NAME
                + ") must be at least one greater than genotype calculation max alternate alleles ("
                + GenotypeCalculationArgumentCollection.MAX_ALTERNATE_ALLELES_LONG_NAME + "), accounting for the non-ref allele");
        }
        this.maxDiploidAltAllelesThatCanBeGenotyped = genomicsdbArgs.maxDiploidAltAllelesThatCanBeGenotyped;
        if (genotypeCalcArgs != null) {
            this.maxGenotypeCount = genotypeCalcArgs.maxGenotypeCount;
        } else {
            // Some defaults
            this.maxGenotypeCount = GenotypeCalculationArgumentCollection.DEFAULT_MAX_GENOTYPE_COUNT;
        }
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
