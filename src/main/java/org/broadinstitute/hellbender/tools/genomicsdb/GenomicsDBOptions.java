package org.broadinstitute.hellbender.tools.genomicsdb;

import java.nio.file.Path;

/**
 * Encapsulates the GenomicsDB-specific options relevant to the FeatureDataSource
 */
public final class GenomicsDBOptions {
    final private Path reference;
    final private boolean callGenotypes;

    public GenomicsDBOptions() {
        this(null, false);
    }

    public GenomicsDBOptions(final Path reference) {
        this(reference, false);
    }

    /**
     *
     * @param reference     Path to a reference. May be null. Needed only for reading from GenomicsDB.
     * @param callGenotypes Indicated whether GenomicsDB should return called genotypes
     */
    public GenomicsDBOptions(final Path reference, final boolean callGenotypes) {
        this.reference = reference;
        this.callGenotypes = callGenotypes;
    }

    public Path getReference() {
        return reference;
    }

    public boolean doCallGenotypes() {
        return callGenotypes;
    }
}
