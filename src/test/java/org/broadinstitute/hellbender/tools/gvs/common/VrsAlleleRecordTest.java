package org.broadinstitute.hellbender.tools.gvs.common;

import org.testng.annotations.Test;

import static org.testng.Assert.*;

/**
 * Unit tests for VrsAlleleRecord.
 * 
 * VrsAlleleRecord is a simple transport object holding just the VRS allele ID.
 * Location information is stored in the canonical vrs_allele table.
 */
public class VrsAlleleRecordTest {
    
    private static final String TEST_ALLELE_ID = "ga4gh:VA.abc123def456";
    
    @Test
    public void testConstructor() {
        VrsAlleleRecord record = new VrsAlleleRecord(TEST_ALLELE_ID);
        assertEquals(record.vrsAlleleId, TEST_ALLELE_ID);
    }
    
    @Test
    public void testMultipleAlleleIds() {
        VrsAlleleRecord rec0 = new VrsAlleleRecord("ga4gh:VA.id0");
        VrsAlleleRecord rec1 = new VrsAlleleRecord("ga4gh:VA.id1");
        VrsAlleleRecord rec2 = new VrsAlleleRecord("ga4gh:VA.id2");
        
        assertEquals(rec0.vrsAlleleId, "ga4gh:VA.id0");
        assertEquals(rec1.vrsAlleleId, "ga4gh:VA.id1");
        assertEquals(rec2.vrsAlleleId, "ga4gh:VA.id2");
    }
    
    @Test
    public void testEquality() {
        VrsAlleleRecord rec1 = new VrsAlleleRecord(TEST_ALLELE_ID);
        VrsAlleleRecord rec2 = new VrsAlleleRecord(TEST_ALLELE_ID);
        
        assertEquals(rec1, rec2);
    }
    
    @Test
    public void testInequalityByAlleleId() {
        VrsAlleleRecord rec1 = new VrsAlleleRecord("ga4gh:VA.id1");
        VrsAlleleRecord rec2 = new VrsAlleleRecord("ga4gh:VA.id2");
        
        assertNotEquals(rec1, rec2);
    }
    
    @Test
    public void testHashCodeConsistency() {
        VrsAlleleRecord rec1 = new VrsAlleleRecord(TEST_ALLELE_ID);
        VrsAlleleRecord rec2 = new VrsAlleleRecord(TEST_ALLELE_ID);
        
        // Equal objects must have equal hash codes
        assertEquals(rec1.hashCode(), rec2.hashCode());
    }
    
    @Test
    public void testHashCodeDifference() {
        VrsAlleleRecord rec1 = new VrsAlleleRecord("ga4gh:VA.id1");
        VrsAlleleRecord rec2 = new VrsAlleleRecord("ga4gh:VA.id2");
        
        // Different allele IDs should have different hash codes
        assertNotEquals(rec1.hashCode(), rec2.hashCode());
    }
    
    @Test
    public void testToString() {
        VrsAlleleRecord record = new VrsAlleleRecord(TEST_ALLELE_ID);
        String str = record.toString();
        
        // Verify key components are in the string representation
        assertTrue(str.contains(TEST_ALLELE_ID));
    }
    
    @Test
    public void testHetVarScenario() {
        // Simulate het-var with genotype 1/2
        // Variant: ref="C", alt="T,ATG" at position 12345
        // Two alternate alleles, each gets its own record
        
        VrsAlleleRecord alt0 = new VrsAlleleRecord("ga4gh:VA.forT");
        VrsAlleleRecord alt1 = new VrsAlleleRecord("ga4gh:VA.forATG");
        
        // Different allele IDs
        assertNotEquals(alt0.vrsAlleleId, alt1.vrsAlleleId);
        
        // Different records (not equal)
        assertNotEquals(alt0, alt1);
    }
    
    @Test
    public void testSelfEquality() {
        VrsAlleleRecord record = new VrsAlleleRecord(TEST_ALLELE_ID);
        assertEquals(record, record);
    }
    
    @Test
    public void testNotEqualToNull() {
        VrsAlleleRecord record = new VrsAlleleRecord(TEST_ALLELE_ID);
        assertNotEquals(record, null);
    }
    
    @Test
    public void testNotEqualToOtherType() {
        VrsAlleleRecord record = new VrsAlleleRecord(TEST_ALLELE_ID);
        assertNotEquals(record, "not a VrsAlleleRecord");
    }
    
    @Test
    public void testMultiAllelicVariant() {
        // Multi-allelic variant: ref="A", alt="T,TT,TTT"
        // Three alternate alleles
        VrsAlleleRecord[] records = new VrsAlleleRecord[]{
            new VrsAlleleRecord("ga4gh:VA.id1"),
            new VrsAlleleRecord("ga4gh:VA.id2"),
            new VrsAlleleRecord("ga4gh:VA.id3")
        };
        
        assertEquals(records.length, 3);
        assertNotEquals(records[0], records[1]);
        assertNotEquals(records[1], records[2]);
        assertNotEquals(records[0], records[2]);
    }
}

