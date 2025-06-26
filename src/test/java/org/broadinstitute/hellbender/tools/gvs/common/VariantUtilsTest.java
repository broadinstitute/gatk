package org.broadinstitute.hellbender.tools.gvs.common;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class VariantUtilsTest {

    @Test
    public void testRightTrimBasicExample1() {
        // Test case 1 from the prompt
        final String reference = "CACGTACGT";
        final String allele = "ACGTACGT";
        final Pair<String, String> expected = Pair.of("CA", "A");

        assertEquals(VariantUtils.rightTrim(reference, allele), expected);
    }

    @Test
    public void testRightTrimBasicExample2() {
        // Test case 2 from the prompt
        final String reference = "ACGTACGT";
        final String allele = "GGACGTACGT";
        final Pair<String, String> expected = Pair.of("A", "GGA");

        assertEquals(VariantUtils.rightTrim(reference, allele), expected);
    }

    @Test
    public void testRightTrimNoCommonSuffix() {
        // No common suffix, sequences should remain unchanged
        final String reference = "AAAA";
        final String allele = "TTTT";
        final Pair<String, String> expected = Pair.of("AAAA", "TTTT");

        assertEquals(VariantUtils.rightTrim(reference, allele), expected);
    }

    @Test
    public void testRightTrimSingleCharacters() {
        // Single character sequences should remain unchanged
        final String reference = "A";
        final String allele = "T";
        final Pair<String, String> expected = Pair.of("A", "T");

        assertEquals(VariantUtils.rightTrim(reference, allele), expected);
    }

    @Test
    public void testRightTrimPartialCommonSuffix() {
        // Partial common suffix
        final String reference = "ATCGATCG";
        final String allele = "GGATCG";
        final Pair<String, String> expected = Pair.of("ATC", "G");

        assertEquals(VariantUtils.rightTrim(reference, allele), expected);
    }

    @Test
    public void testRightTrimDifferentLengths() {
        // Different length sequences with common suffix
        final String reference = "AAATTTCCC";
        final String allele = "GCCC";
        final Pair<String, String> expected = Pair.of("AAATTT", "G");

        assertEquals(VariantUtils.rightTrim(reference, allele), expected);
    }

    @Test
    public void testRightTrimLongCommonSuffix() {
        // Long common suffix that gets truncated to preserve minimum length
        final String reference = "AATCGATCGATCG";
        final String allele = "GGTCGATCGATCG";
        final Pair<String, String> expected = Pair.of("AA", "GG");

        assertEquals(VariantUtils.rightTrim(reference, allele), expected);
    }

    @Test
    public void testRightTrimOneCharacterWithSuffix() {
        // One sequence is single character, other has matching suffix
        final String reference = "G";
        final String allele = "TTG";
        final Pair<String, String> expected = Pair.of("G", "TTG");

        assertEquals(VariantUtils.rightTrim(reference, allele), expected);
    }

    @Test
    public void testRightTrimInsertion() {
        // Simulates an insertion variant
        final String reference = "A";
        final String allele = "ATCG";
        final Pair<String, String> expected = Pair.of("A", "ATCG");

        assertEquals(VariantUtils.rightTrim(reference, allele), expected);
    }

    @Test
    public void testRightTrimDeletion() {
        // Simulates a deletion variant
        final String reference = "ATCG";
        final String allele = "A";
        final Pair<String, String> expected = Pair.of("ATCG", "A");

        assertEquals(VariantUtils.rightTrim(reference, allele), expected);
    }

    @Test
    public void testRightTrimComplexVariant() {
        // Complex variant with substitution and common suffix
        final String reference = "ATCGATCGATCGATCG";
        final String allele = "GGCGATCGATCGATCG";
        final Pair<String, String> expected = Pair.of("AT", "GG");

        assertEquals(VariantUtils.rightTrim(reference, allele), expected);
    }

    @Test(expectedExceptions = UserException.class)
    public void testRightTrimNullReference() {
        VariantUtils.rightTrim(null, "ATCG");
    }

    @Test(expectedExceptions = UserException.class)
    public void testRightTrimNullAllele() {
        VariantUtils.rightTrim("ATCG", null);
    }

    @Test(expectedExceptions = UserException.class)
    public void testRightTrimEmptyReference() {
        VariantUtils.rightTrim("", "ATCG");
    }

    @Test(expectedExceptions = UserException.class)
    public void testRightTrimEmptyAllele() {
        VariantUtils.rightTrim("ATCG", "");
    }

    @Test(expectedExceptions = UserException.class)
    public void testRightTrimBothEmpty() {
        VariantUtils.rightTrim("", "");
    }

    @Test(expectedExceptions = UserException.class)
    public void testRightTrimBothNull() {
        VariantUtils.rightTrim(null, null);
    }
}
