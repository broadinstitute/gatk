package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * Test class for {@link DataSourceUtils}
 * Created by jonn on 3/22/18.
 */
public class DataSourceUtilsUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Iterator<Object[]> provideForValidateVersionInformation() {
        final ArrayList<Object[]> testArgs = new ArrayList<>();

        final ArrayList<Integer> baseArgs = new ArrayList<>();
        baseArgs.add(DataSourceUtils.MIN_MAJOR_VERSION_NUMBER);
        baseArgs.add(DataSourceUtils.MIN_MINOR_VERSION_NUMBER);
        baseArgs.add(DataSourceUtils.MIN_YEAR_RELEASED);
        baseArgs.add(DataSourceUtils.MIN_MONTH_RELEASED);
        baseArgs.add(DataSourceUtils.MIN_DAY_RELEASED);

        for ( int offset = -1 ; offset < 2; ++offset ) {
            for ( int i = 0; i < baseArgs.size(); ++i ) {

                final ArrayList<Object> argList = new ArrayList<>();

                argList.addAll(baseArgs.subList(0, i));
                argList.add(baseArgs.get(i) + offset);

                if ( i < baseArgs.size() - 1) {
                    argList.addAll(baseArgs.subList(i + 1, baseArgs.size()));
                }

                argList.add(offset >= 0);

                testArgs.add(argList.toArray());
            }
        }

        return testArgs.iterator();
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForValidateVersionInformation")
    public void testValidateVersionInformation(final Integer major,
                                               final Integer minor,
                                               final Integer year,
                                               final Integer month,
                                               final Integer day,
                                               final Boolean expected) {
        Assert.assertEquals(
                DataSourceUtils.validateVersionInformation(
                    major.intValue(),
                    minor.intValue(),
                    year.intValue(),
                    month.intValue(),
                    day.intValue()),
            expected.booleanValue() );
    }

}
