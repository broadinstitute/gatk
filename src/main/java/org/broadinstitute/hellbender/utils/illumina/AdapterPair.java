package org.broadinstitute.hellbender.utils.illumina;

public interface AdapterPair {
    String get3PrimeAdapter();
    String get5PrimeAdapter();
    String getName();
}
