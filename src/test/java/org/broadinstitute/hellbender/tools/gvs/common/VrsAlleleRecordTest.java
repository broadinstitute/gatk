package org.broadinstitute.hellbender.tools.gvs.common;

import org.testng.annotations.Test;

import static org.testng.Assert.*;

/**
 * Unit tests for VrsAlleleRecord.
 *
 * VrsAlleleRecord is a full canonical VRS record carrying all columns needed
 * to populate the vrs_allele table.
 */
public class VrsAlleleRecordTest {

    // Shared test constants
    private static final String ALLELE_ID   = "ga4gh:VA.abc123def456";
    private static final String LOCATION_ID = "ga4gh:SL.loc001";
    private static final String REFGET      = "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl";
    private static final long   START       = 12344L;  // 0-based interbase
    private static final long   END         = 12345L;
    private static final String ALT_SEQ     = "T";

    // ------------------------------------------------------------------
    // LiteralSequenceExpression constructor
    // ------------------------------------------------------------------

    @Test
    public void testLiteralConstructorFields() {
        VrsAlleleRecord rec = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, START, END, ALT_SEQ);
        assertEquals(rec.vrsAlleleId,    ALLELE_ID);
        assertEquals(rec.vrsLocationId,  LOCATION_ID);
        assertEquals(rec.refgetAccession, REFGET);
        assertEquals(rec.start, START);
        assertEquals(rec.end,   END);
        assertEquals(rec.stateType,     VrsAlleleRecord.STATE_LITERAL);
        assertEquals(rec.stateSequence, ALT_SEQ);
        assertNull(rec.stateLength);
        assertNull(rec.stateRepeatSubunitLength);
    }

    @Test
    public void testLiteralStateTypeConstant() {
        VrsAlleleRecord rec = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, START, END, ALT_SEQ);
        assertEquals(rec.stateType, "LiteralSequenceExpression");
    }

    // ------------------------------------------------------------------
    // ReferenceLengthExpression constructor
    // ------------------------------------------------------------------

    @Test
    public void testRleConstructorFields() {
        VrsAlleleRecord rec = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, 100L, 106L, 3, 3);
        assertEquals(rec.vrsAlleleId,    ALLELE_ID);
        assertEquals(rec.vrsLocationId,  LOCATION_ID);
        assertEquals(rec.refgetAccession, REFGET);
        assertEquals(rec.start, 100L);
        assertEquals(rec.end,   106L);
        assertEquals(rec.stateType, VrsAlleleRecord.STATE_REF_LEN);
        assertNull(rec.stateSequence);
        assertEquals(rec.stateLength,             Integer.valueOf(3));
        assertEquals(rec.stateRepeatSubunitLength, Integer.valueOf(3));
    }

    @Test
    public void testRleStateTypeConstant() {
        VrsAlleleRecord rec = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, 100L, 106L, 3, 3);
        assertEquals(rec.stateType, "ReferenceLengthExpression");
    }

    @Test
    public void testRleWithDifferentLengthAndSubunit() {
        // deletion: length=0, repeatSubunitLength=6
        VrsAlleleRecord rec = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, 200L, 206L, 0, 6);
        assertEquals(rec.stateLength,             Integer.valueOf(0));
        assertEquals(rec.stateRepeatSubunitLength, Integer.valueOf(6));
    }

    // ------------------------------------------------------------------
    // Equality and hash code
    // ------------------------------------------------------------------

    @Test
    public void testLiteralEquality() {
        VrsAlleleRecord rec1 = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, START, END, ALT_SEQ);
        VrsAlleleRecord rec2 = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, START, END, ALT_SEQ);
        assertEquals(rec1, rec2);
        assertEquals(rec1.hashCode(), rec2.hashCode());
    }

    @Test
    public void testRleEquality() {
        VrsAlleleRecord rec1 = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, 100L, 106L, 3, 3);
        VrsAlleleRecord rec2 = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, 100L, 106L, 3, 3);
        assertEquals(rec1, rec2);
        assertEquals(rec1.hashCode(), rec2.hashCode());
    }

    @Test
    public void testInequalityByAlleleId() {
        VrsAlleleRecord rec1 = new VrsAlleleRecord("ga4gh:VA.id1", LOCATION_ID, REFGET, START, END, ALT_SEQ);
        VrsAlleleRecord rec2 = new VrsAlleleRecord("ga4gh:VA.id2", LOCATION_ID, REFGET, START, END, ALT_SEQ);
        assertNotEquals(rec1, rec2);
    }

    @Test
    public void testInequalityByStateType() {
        VrsAlleleRecord literal = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, START, END, ALT_SEQ);
        VrsAlleleRecord rle     = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, START, END, 1, 1);
        assertNotEquals(literal, rle);
    }

    @Test
    public void testSelfEquality() {
        VrsAlleleRecord rec = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, START, END, ALT_SEQ);
        assertEquals(rec, rec);
    }

    @Test
    public void testNotEqualToNull() {
        VrsAlleleRecord rec = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, START, END, ALT_SEQ);
        assertNotEquals(rec, null);
    }

    @Test
    public void testNotEqualToOtherType() {
        VrsAlleleRecord rec = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, START, END, ALT_SEQ);
        assertNotEquals(rec, "not a VrsAlleleRecord");
    }

    // ------------------------------------------------------------------
    // toString
    // ------------------------------------------------------------------

    @Test
    public void testToStringContainsKeyFields() {
        VrsAlleleRecord rec = new VrsAlleleRecord(ALLELE_ID, LOCATION_ID, REFGET, START, END, ALT_SEQ);
        String str = rec.toString();
        assertTrue(str.contains(ALLELE_ID));
        assertTrue(str.contains(LOCATION_ID));
        assertTrue(str.contains(REFGET));
    }

    // ------------------------------------------------------------------
    // Multi-allelic scenario
    // ------------------------------------------------------------------

    @Test
    public void testHetVarScenario() {
        // Simulate het-var with genotype 1/2
        // Variant: ref="C", alt="T,ATG" at position 12345
        VrsAlleleRecord alt0 = new VrsAlleleRecord("ga4gh:VA.forT",   LOCATION_ID, REFGET, START, END, "T");
        VrsAlleleRecord alt1 = new VrsAlleleRecord("ga4gh:VA.forATG", LOCATION_ID, REFGET, START, 12345L, "ATG");

        assertNotEquals(alt0.vrsAlleleId, alt1.vrsAlleleId);
        assertNotEquals(alt0, alt1);
    }

    @Test
    public void testMultiAllelicVariant() {
        // Multi-allelic variant: ref="A", alt="T,TT,TTT"
        VrsAlleleRecord[] records = {
            new VrsAlleleRecord("ga4gh:VA.id1", LOCATION_ID, REFGET, START, END, "T"),
            new VrsAlleleRecord("ga4gh:VA.id2", LOCATION_ID, REFGET, START, END, "TT"),
            new VrsAlleleRecord("ga4gh:VA.id3", LOCATION_ID, REFGET, START, END, "TTT"),
        };

        assertEquals(records.length, 3);
        assertNotEquals(records[0], records[1]);
        assertNotEquals(records[1], records[2]);
        assertNotEquals(records[0], records[2]);
    }
}
