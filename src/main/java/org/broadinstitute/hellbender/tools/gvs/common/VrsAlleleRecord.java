package org.broadinstitute.hellbender.tools.gvs.common;

/**
 * Represents a normalized VRS allele with its computed identifier.
 * 
 * Simple transport object used during VET ingest to collect VRS allele IDs
 * for each alternate allele in a variant context. The IDs are then assembled
 * into an array for storage in the VET table.
 * 
 * Location information is NOT stored here; it is part of the vrs_allele
 * table in the canonical VRS schema.
 */
public class VrsAlleleRecord {
    
    /**
     * VRS Allele identifier (ga4gh:VA.<digest>).
     * Computed from the normalized (start, end, altSequence) representation.
     * Location information is stored separately in the canonical vrs_allele table.
     */
    public final String vrsAlleleId;
    
    /**
     * Constructor.
     *
     * @param vrsAlleleId "ga4gh:VA.<digest>"
     */
    public VrsAlleleRecord(String vrsAlleleId) {
        this.vrsAlleleId = vrsAlleleId;
    }
    
    @Override
    public String toString() {
        return String.format("VrsAlleleRecord{vrsAlleleId='%s'}", vrsAlleleId);
    }
    
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof VrsAlleleRecord)) return false;
        
        VrsAlleleRecord that = (VrsAlleleRecord) o;
        return vrsAlleleId.equals(that.vrsAlleleId);
    }
    
    @Override
    public int hashCode() {
        return vrsAlleleId.hashCode();
    }
}
