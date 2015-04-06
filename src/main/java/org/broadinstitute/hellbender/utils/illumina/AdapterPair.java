package org.broadinstitute.hellbender.utils.illumina;

public interface AdapterPair {

    public String get3PrimeAdapter();
    public String get3PrimeAdapterInReadOrder();
    public byte[] get3PrimeAdapterBytes();
    public byte[] get3PrimeAdapterBytesInReadOrder();

    public String get5PrimeAdapter();
    public String get5PrimeAdapterInReadOrder();
    public byte[] get5PrimeAdapterBytes();
    public byte[] get5PrimeAdapterBytesInReadOrder();

    public String getName();
}
