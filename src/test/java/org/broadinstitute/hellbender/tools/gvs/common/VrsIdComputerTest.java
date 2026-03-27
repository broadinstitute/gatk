package org.broadinstitute.hellbender.tools.gvs.common;

import org.testng.annotations.Test;

import static org.testng.Assert.assertEquals;

/**
 * Unit tests for VrsIdComputer.
 * 
 * Tests cover:
 * - SequenceLocation ID computation
 * - Allele ID computation for various variant types
 * - Digest stability (canonical JSON output)
 * - Format validation (ga4gh: prefix, 32-char digest)
 */
public class VrsIdComputerTest {

    private static final String REFGET_SQ = "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl";

    @Test
    public void testLocationIdFormat() {
        // Verify output format: ga4gh:SL.{32-char digest}
        String locId = VrsIdComputer.computeLocationId(REFGET_SQ, 0, 1);
        assertEquals(locId.substring(0, 9), "ga4gh:SL.");
        assertEquals(locId.length(), 41); // 9 prefix + 32 digest = 41
    }

    @Test
    public void testAlleleIdFormat() {
        // Verify output format: ga4gh:VA.{32-char digest}
        String alleleId = VrsIdComputer.computeAlleleId(REFGET_SQ, 0, 1, "T");
        assertEquals(alleleId.substring(0, 9), "ga4gh:VA.");
        assertEquals(alleleId.length(), 41); // 9 prefix + 32 digest = 41
    }

    @Test
    public void testLocationIdDeterministic() {
        // Same inputs should produce same ID (deterministic JSON serialization)
        String locId1 = VrsIdComputer.computeLocationId(REFGET_SQ, 10, 15);
        String locId2 = VrsIdComputer.computeLocationId(REFGET_SQ, 10, 15);
        assertEquals(locId1, locId2);
    }

    @Test
    public void testAlleleIdDeterministic() {
        // Same inputs should produce same ID
        String id1 = VrsIdComputer.computeAlleleId(REFGET_SQ, 10, 15, "ATG");
        String id2 = VrsIdComputer.computeAlleleId(REFGET_SQ, 10, 15, "ATG");
        assertEquals(id1, id2);
    }

    @Test
    public void testLocationIdVariesByCoordinates() {
        // Different coordinates should produce different IDs
        String loc1 = VrsIdComputer.computeLocationId(REFGET_SQ, 0, 1);
        String loc2 = VrsIdComputer.computeLocationId(REFGET_SQ, 0, 2);
        String loc3 = VrsIdComputer.computeLocationId(REFGET_SQ, 1, 2);
        
        // All three should be different
        assertEquals(loc1.equals(loc2), false);
        assertEquals(loc1.equals(loc3), false);
        assertEquals(loc2.equals(loc3), false);
    }

    @Test
    public void testAlleleIdVariesByAltSequence() {
        // Different alt sequences should produce different IDs (same location)
        String id1 = VrsIdComputer.computeAlleleId(REFGET_SQ, 10, 11, "A");
        String id2 = VrsIdComputer.computeAlleleId(REFGET_SQ, 10, 11, "T");
        String id3 = VrsIdComputer.computeAlleleId(REFGET_SQ, 10, 11, "ATG");
        
        // All three should be different
        assertEquals(id1.equals(id2), false);
        assertEquals(id1.equals(id3), false);
        assertEquals(id2.equals(id3), false);
    }

    @Test
    public void testAlleleIdVariesByLocation() {
        // Same alt sequence, different locations should produce different IDs
        String id1 = VrsIdComputer.computeAlleleId(REFGET_SQ, 0, 1, "T");
        String id2 = VrsIdComputer.computeAlleleId(REFGET_SQ, 0, 2, "T");
        String id3 = VrsIdComputer.computeAlleleId(REFGET_SQ, 1, 2, "T");
        
        // All three should be different
        assertEquals(id1.equals(id2), false);
        assertEquals(id1.equals(id3), false);
        assertEquals(id2.equals(id3), false);
    }

    @Test
    public void testAlleleIdVariesByRefgetAccession() {
        // Different RefGet accession should produce different ID
        String id1 = VrsIdComputer.computeAlleleId(REFGET_SQ, 10, 11, "A");
        String id2 = VrsIdComputer.computeAlleleId("SQ.DifferentAccession", 10, 11, "A");
        
        assertEquals(id1.equals(id2), false);
    }

    @Test
    public void testSnpAlleleId() {
        // SNP: substitution from ref→alt
        // Position 10, single nucleotide substitution (start=10, end=11)
        String snpId = VrsIdComputer.computeAlleleId(REFGET_SQ, 10, 11, "G");
        assertEquals(snpId.substring(0, 9), "ga4gh:VA.");
        assertEquals(snpId.length(), 41);
    }

    @Test
    public void testDeletionAlleleId() {
        // Deletion: alt sequence is empty string
        // Deletion of 5 bp at position 20-25
        String delId = VrsIdComputer.computeAlleleId(REFGET_SQ, 20, 25, "");
        assertEquals(delId.substring(0, 9), "ga4gh:VA.");
        assertEquals(delId.length(), 41);
    }

    @Test
    public void testInsertionAlleleId() {
        // Insertion: alt sequence inserted at position (ref is empty)
        // Insertion of "ATGC" at position 50
        String insId = VrsIdComputer.computeAlleleId(REFGET_SQ, 50, 50, "ATGC");
        assertEquals(insId.substring(0, 9), "ga4gh:VA.");
        assertEquals(insId.length(), 41);
    }

    @Test
    public void testComplexAlleleId() {
        // Complex variant: multi-nucleotide change
        String complexId = VrsIdComputer.computeAlleleId(REFGET_SQ, 100, 105, "ATCG");
        assertEquals(complexId.substring(0, 9), "ga4gh:VA.");
        assertEquals(complexId.length(), 41);
    }

    @Test
    public void testLocationIdZeroBased() {
        // Test with 0-based coordinates (VRS standard)
        String loc0 = VrsIdComputer.computeLocationId(REFGET_SQ, 0, 0);
        String loc1 = VrsIdComputer.computeLocationId(REFGET_SQ, 0, 1);
        
        // Empty range vs. 1 bp range should differ
        assertEquals(loc0.equals(loc1), false);
    }

    @Test
    public void testLocationIdLargeCoordinates() {
        // Test with large genomic coordinates (chromosome-scale)
        String locLarge = VrsIdComputer.computeLocationId(REFGET_SQ, 1_000_000, 1_000_100);
        assertEquals(locLarge.substring(0, 9), "ga4gh:SL.");
        assertEquals(locLarge.length(), 41);
    }

    @Test
    public void testAlleleIdSpecialCharacters() {
        // Test alt sequence with all valid DNA bases
        String idDna = VrsIdComputer.computeAlleleId(REFGET_SQ, 0, 1, "ATCGN");
        assertEquals(idDna.substring(0, 9), "ga4gh:VA.");
        assertEquals(idDna.length(), 41);
    }

    @Test
    public void testRefgetAccessionWithSpecialCharacters() {
        // RefGet accessions can contain '_' and numbers
        String id = VrsIdComputer.computeAlleleId("SQ.abc123_xyz", 10, 11, "A");
        assertEquals(id.substring(0, 9), "ga4gh:VA.");
        assertEquals(id.length(), 41);
    }

    @Test
    public void testAlleleAndLocationConsistency() {
        // Allele ID computation should include LocationId, so changing location changes Allele ID
        String allele1 = VrsIdComputer.computeAlleleId(REFGET_SQ, 10, 11, "A");
        String allele2 = VrsIdComputer.computeAlleleId(REFGET_SQ, 10, 12, "A");
        
        // Even with same alt, different locations produce different Allele IDs
        assertEquals(allele1.equals(allele2), false);
    }

    // ── ReferenceLengthExpression tests ───────────────────────────────────────

    @Test
    public void testReferenceLengthAlleleIdFormat() {
        // Verify output format: ga4gh:VA.{32-char digest}
        String id = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(REFGET_SQ, 20, 25, 5, 5);
        assertEquals(id.substring(0, 9), "ga4gh:VA.");
        assertEquals(id.length(), 41);
    }

    @Test
    public void testReferenceLengthAlleleIdDeterministic() {
        // Same inputs must produce the same ID every time
        String id1 = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(REFGET_SQ, 20, 25, 5, 5);
        String id2 = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(REFGET_SQ, 20, 25, 5, 5);
        assertEquals(id1, id2);
    }

    @Test
    public void testReferenceLengthDiffersFromLiteralForSameLocation() {
        // A deletion represented as ReferenceLengthExpression must NOT equal the same
        // location with a LiteralSequenceExpression (different state type → different ID)
        String refLengthId = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(REFGET_SQ, 20, 25, 5, 5);
        String literalId   = VrsIdComputer.computeAlleleId(REFGET_SQ, 20, 25, "ATCGA");
        assertEquals(refLengthId.equals(literalId), false);
    }

    @Test
    public void testReferenceLengthVariesByLength() {
        // Deletions of different sizes must produce different IDs
        String id5 = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(REFGET_SQ, 20, 25, 5, 5);
        String id3 = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(REFGET_SQ, 20, 23, 3, 3);
        assertEquals(id5.equals(id3), false);
    }

    @Test
    public void testReferenceLengthVariesByRepeatSubunitLength() {
        // Same total length but different repeat unit length (deletion vs tandem repeat)
        // length=6, repeatSubunitLength=6 (simple deletion) vs repeatSubunitLength=2 (repeat unit)
        String simpleDel  = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(REFGET_SQ, 20, 26, 6, 6);
        String repeatDel  = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(REFGET_SQ, 20, 26, 6, 2);
        assertEquals(simpleDel.equals(repeatDel), false);
    }

    @Test
    public void testReferenceLengthVariesByLocation() {
        // Same deletion size at different positions
        String id1 = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(REFGET_SQ, 20, 25, 5, 5);
        String id2 = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(REFGET_SQ, 30, 35, 5, 5);
        assertEquals(id1.equals(id2), false);
    }

    @Test
    public void testReferenceLengthVariesByRefgetAccession() {
        // Same deletion, different chromosome → different ID
        String id1 = VrsIdComputer.computeAlleleIdForReferenceLengthExpression(REFGET_SQ, 20, 25, 5, 5);
        String id2 = VrsIdComputer.computeAlleleIdForReferenceLengthExpression("SQ.DifferentAccession", 20, 25, 5, 5);
        assertEquals(id1.equals(id2), false);
    }
}
