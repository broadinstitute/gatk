package org.broadinstitute.hellbender.utils;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public final class FisherExactTestUnitTest{
    private static final double DELTA_PRECISION = 10e-7;

    @DataProvider(name = "UsingTable")
    public Object[][] makeUsingTable() {

        /* NOTE: the expected P values were computed in R as follows
             fisher = function(v) {
                return(fisher.test(matrix(v, nrow=2, ncol=2))$p.value)
             }
        */
        //> fisher(c(2068, 6796, 1133, 0))

        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{0, 0, 0, 0,                     1.0});
        tests.add(new Object[]{100000, 100000, 100000, 100000, 1.0}  );
        tests.add(new Object[]{1, 2, 3, 4,                     1.0});
        tests.add(new Object[]{0, 0, 100000, 100000,           1.0});
        tests.add(new Object[]{100000,100000,100000,0,         0.0});             //below R's or Java's precision
        tests.add(new Object[]{12, 4, 26, 7,                   1.0});             //tests rounding the probabilities from the Hypergeometric
        tests.add(new Object[]{12, 26, 4, 7,                   1.0});             //tests rounding the probabilities from the Hypergeometric

        tests.add(new Object[]{200000,100000,1,2,              0.259263});
        tests.add(new Object[]{100,100,100,0,                  3.730187e-23});
        tests.add(new Object[]{13736,9047,41,1433,             0.0});              //below R's or Java's precision
        tests.add(new Object[]{66, 14, 64, 4,                  0.0424333});
        tests.add(new Object[]{351169, 306836, 153739, 2379,   0.0});              //below R's or Java's precision
        tests.add(new Object[]{116449, 131216, 289, 16957,     0.0});              //below R's or Java's precision
        tests.add(new Object[]{137, 159, 9, 23,                0.06088506});
        tests.add(new Object[]{129, 90, 21, 20,                0.3919603});
        tests.add(new Object[]{14054, 9160, 16, 7827,          0.0});             //below R's or Java's precision
        tests.add(new Object[]{32803, 9184, 32117, 3283,       0.0});             //below R's or Java's precision
        tests.add(new Object[]{2068, 6796, 1133, 0,            0.0});             //below R's or Java's precision

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "UsingTable")
    public void testUsingTableData(final int refpos, final int refneg, final int altpos, final int altneg, double expectedPvalue) {
        final int[][] contingencyTable = new int[2][2];
        contingencyTable[0][0] = refpos;
        contingencyTable[0][1] = refneg;
        contingencyTable[1][0] = altpos;
        contingencyTable[1][1] = altneg;
        final double pvalue = FisherExactTest.twoSidedPValue(contingencyTable);
        Assert.assertEquals(pvalue, expectedPvalue, DELTA_PRECISION, "Pvalues");
    }
}
