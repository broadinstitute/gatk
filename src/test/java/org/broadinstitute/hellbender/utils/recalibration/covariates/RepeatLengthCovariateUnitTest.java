package org.broadinstitute.hellbender.utils.recalibration.covariates;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class RepeatLengthCovariateUnitTest extends GATKBaseTest {

    @DataProvider
    public static Object[][] FindTandemRepeatUnitsTestData() {
        return new Object[][]{
                // Parentheses added for readability (removed in the test).
                // Current implementation: the repeat unit is determined such that the character at the offset
                // will be the *last* character in the repeat unit in the case of two letter repeats.
                // In other words the repeat unit is determined up to cyclic permutations.

                // { readBaseString, offset, expectedRepeatUnit, expectedRepeatLength }
                {"A_CA_C(A)_CA_CA_C", 4, "CA", 4},
                {"(A)_CA_CA_CA_CA_C", 0, "CA", 4},
                {"(A)_CTA_CTA_CTA_CTA_CT", 0, "CTA", 4},
                {"A_CTA_CT(A)_CTA_CTA_CT", 6, "CTA", 4},
                {"AC_TAC_TA(C)_TAC_TAC_T", 7, "TAC", 4},
                {"ACT_ACT_AC(T)_ACT_ACT", 8, "ACT", 5},
                {"A_CTA_GCT_AC(T)_ACT_ACT", 9, "ACT", 3},
                {"(A)AAAAAAAAA", 0, "A", 10},
                {"AAAA(A)AAAAA", 4, "A", 10},
                {"A_C_(A)AAAAA", 2, "A", 6},
                {"AAAAAAAAA(A)", 9, "A", 10},
                {"A(C)_TACT_TACT_ACTA_CTAC", 1, "TACT", 2},
                // this (ACT)*2 is missed. When max repeat count is 1, the base at offset+1 is chosen as the repeat unit.
                {"AC(T)_ACT_TAC_TAC_TAC_TAC", 2, "A", 1},
                {"TTTTTT(C)_AAA_AAA_AAA", 6, "A", 9}, // Note that the repeat unit is "A", not "AAA"
                {"TTT(T)TTC_AAA_AAA_AAA", 3, "T", 6},
                {"A_CTA_CT(G)_ACT_TAC_TAC_TAC_TAC", 6, "A", 1}, // "TAC" downstream is not recognized.
                {"CT_ACT_ACT_A(C)T_ACT_ACT", 9, "TAC", 5}, // "ACT" is not detected, but "TAC" is, which is the same repeat up to cyclic permutation.
                {"ATTT_(A)TTT_ATTT_CTTT", 4, "T", 3}
        };
    }

    @Test(dataProvider="FindTandemRepeatUnitsTestData")
    public void testFindTandemRepeatUnits(final String readBaseString, final int offset,
                                          final String expectedRepeatUnit, final int expectedRepeatLength){
        // remove the characters that were put in to make the input string easier to read
        final byte[] readBases = readBaseString.replaceAll("[()_]", "").getBytes();
        RepeatLengthCovariate repeatLengthCovariate = new RepeatLengthCovariate();
        // tsato: do I really need to initialize this...
        repeatLengthCovariate.initialize(new RecalibrationArgumentCollection(), Arrays.asList("yo"));

        Pair<byte[], Integer> ans = repeatLengthCovariate.findTandemRepeatUnits(readBases, offset);
        // Pair<byte[], Integer> ans = RepeatLengthCovariate.detectSTR(readBases, offset, 8);
        byte[] repeatUnit =  ans.getLeft();
        int repeatLength = ans.getRight();

        // for debugging
        String repeatUnitStr = new String(repeatUnit, StandardCharsets.UTF_8);

        int d = 3;
        Assert.assertEquals(repeatUnitStr, expectedRepeatUnit);
        Assert.assertEquals(repeatLength, expectedRepeatLength);
    }

}
