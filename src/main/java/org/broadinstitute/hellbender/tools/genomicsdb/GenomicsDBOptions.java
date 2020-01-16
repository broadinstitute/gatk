package org.broadinstitute.hellbender.tools.genomicsdb;

import htsjdk.variant.variantcontext.GenotypeLikelihoods;

import java.nio.file.Path;

/**
 * Encapsulates the GenomicsDB-specific options relevant to the FeatureDataSource
 */
public final class GenomicsDBOptions {
    final private Path reference;
    final private boolean callGenotypes;
    final private int maxAlternateAlleles;
    final private int maxGenotypeCount;
    final private int MAX_GENOTYPES = 100000;

    public GenomicsDBOptions() {
        this(null, false, GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);
    }

    public GenomicsDBOptions(final Path reference) {
        this(reference, false, GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);
    }

    public GenomicsDBOptions(final Path reference, final boolean callGenotypes, final int maxAlternateAlleles) {
        this.reference = reference;
        this.callGenotypes = callGenotypes;
        this.maxAlternateAlleles = maxAlternateAlleles;
        this.maxGenotypeCount = MAX_GENOTYPES;
    }

    /**
     *
     * @param reference     Path to a reference. May be null. Needed only for reading from GenomicsDB.
     * @param callGenotypes Indicates whether GenomicsDB should return called genotypes
     * @param maxAlternateAlleles The maximum number of alternative alleles to return in query results
     * @param maxGenotypeCount The maximum number of genotypes to return in query results
     */
    public GenomicsDBOptions(final Path reference, final boolean callGenotypes, final int maxAlternateAlleles, final int maxGenotypeCount) {
        this.reference = reference;
        this.callGenotypes = callGenotypes;
        this.maxAlternateAlleles = maxAlternateAlleles;
        this.maxGenotypeCount = maxGenotypeCount;
    }

    public Path getReference() {
        return reference;
    }

    public boolean doCallGenotypes() {
        return callGenotypes;
    }

    public int getMaxAlternateAlleles() { return maxAlternateAlleles; }

    public int getMaxGenotypeCount() { return maxGenotypeCount; }
}
