package org.broadinstitute.hellbender.tools.gvs.common;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class VariantMinimizationUtilsTest {

    @Test
    public void testMinimizeBasicExample1() {
        // Test case 1 from the prompt
        final String reference = "CACGTACGT";
        final String allele = "ACGTACGT";
        final String expected = "CA,A";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizeBasicExample2() {
        // Test case 2 from the prompt
        final String reference = "ACGTACGT";
        final String allele = "GGACGTACGT";
        final String expected = "A,GGA";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizeNoCommonSuffix() {
        // No common suffix, sequences should remain unchanged
        final String reference = "AAAA";
        final String allele = "TTTT";
        final String expected = "AAAA,TTTT";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizeCompletelyIdentical() {
        // Identical sequences should be minimized to single characters
        final String reference = "ATCG";
        final String allele = "ATCG";
        final String expected = "A,A";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizeSingleCharacters() {
        // Single character sequences should remain unchanged
        final String reference = "A";
        final String allele = "T";
        final String expected = "A,T";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizeSingleCharacterIdentical() {
        // Single identical characters should remain unchanged (minimum length requirement)
        final String reference = "A";
        final String allele = "A";
        final String expected = "A,A";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizePartialCommonSuffix() {
        // Partial common suffix
        final String reference = "ATCGATCG";
        final String allele = "GGATCG";
        final String expected = "ATC,G";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizeDifferentLengths() {
        // Different length sequences with common suffix
        final String reference = "AAATTTCCC";
        final String allele = "GCCC";
        final String expected = "AAATTT,G";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizeLongCommonSuffix() {
        // Long common suffix that gets truncated to preserve minimum length
        final String reference = "AATCGATCGATCG";
        final String allele = "GGTCGATCGATCG";
        final String expected = "AA,GG";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizeOneCharacterWithSuffix() {
        // One sequence is single character, other has matching suffix
        final String reference = "G";
        final String allele = "TTG";
        final String expected = "G,TTG";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizeInsertion() {
        // Simulates an insertion variant
        final String reference = "A";
        final String allele = "ATCG";
        final String expected = "A,ATCG";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizeDeletion() {
        // Simulates a deletion variant
        final String reference = "ATCG";
        final String allele = "A";
        final String expected = "ATCG,A";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test
    public void testMinimizeComplexVariant() {
        // Complex variant with substitution and common suffix
        final String reference = "ATCGATCGATCGATCG";
        final String allele = "GGCGATCGATCGATCG";
        final String expected = "AT,GG";

        assertEquals(VariantMinimizationUtils.minimize(reference, allele), expected);
    }

    @Test(expectedExceptions = UserException.class)
    public void testMinimizeNullReference() {
        VariantMinimizationUtils.minimize(null, "ATCG");
    }

    @Test(expectedExceptions = UserException.class)
    public void testMinimizeNullAllele() {
        VariantMinimizationUtils.minimize("ATCG", null);
    }

    @Test(expectedExceptions = UserException.class)
    public void testMinimizeEmptyReference() {
        VariantMinimizationUtils.minimize("", "ATCG");
    }

    @Test(expectedExceptions = UserException.class)
    public void testMinimizeEmptyAllele() {
        VariantMinimizationUtils.minimize("ATCG", "");
    }

    @Test(expectedExceptions = UserException.class)
    public void testMinimizeBothEmpty() {
        VariantMinimizationUtils.minimize("", "");
    }

    @Test(expectedExceptions = UserException.class)
    public void testMinimizeBothNull() {
        VariantMinimizationUtils.minimize(null, null);
    }
}
