package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import org.apache.commons.lang.RandomStringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;

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

    private List<Object[]>  createBaseTestVersionData() {

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

        // Some specific test cases to prevent regression of version checks:
        // 1 Month after OK release date, but 1 day before OK release date (should pass):
        testArgs.add( new Object[] { DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_YEAR_RELEASED, DataSourceUtils.MIN_MONTH_RELEASED+1, DataSourceUtils.MIN_DAY_RELEASED-1, true } );
        // 1 Year after OK release date, but 1 day before OK release date (should pass):
        testArgs.add( new Object[] { DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_YEAR_RELEASED+1, DataSourceUtils.MIN_MONTH_RELEASED, DataSourceUtils.MIN_DAY_RELEASED-1, true } );
        // 1 Year after OK release date, but 1 month before OK release date (should pass):
        testArgs.add( new Object[] { DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_YEAR_RELEASED+1, DataSourceUtils.MIN_MONTH_RELEASED-1, DataSourceUtils.MIN_DAY_RELEASED, true } );
        // 1 Year after OK release date, but 1 month and 1 day before OK release date (should pass):
        testArgs.add( new Object[] { DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, DataSourceUtils.MIN_YEAR_RELEASED+1, DataSourceUtils.MIN_MONTH_RELEASED-1, DataSourceUtils.MIN_DAY_RELEASED-1, true } );

        return testArgs;
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Iterator<Object[]> provideForValidateVersionInformation() {
        return createBaseTestVersionData().iterator();
    }

    @DataProvider
    private Iterator<Object[]> provideForTestVersionRegex() {

        final ArrayList<Object[]> testArgs = new ArrayList<>();

        final List<Object[]> baseArgs = createBaseTestVersionData();

        for ( final Object[] args : baseArgs ) {
            for ( int whitespace = 0; whitespace < 2; ++whitespace ) {
                for ( int decoratorCount = 0; decoratorCount < 10; ++decoratorCount ) {

                    final String whitespaceString = whitespace != 0 ? "\t \t \t " : " ";
                    final String decoratorString  = decoratorCount  !=0 ? RandomStringUtils.randomAlphanumeric(decoratorCount) : "";

                    testArgs.add(
                            new Object[] {
                                    args[0],args[1],args[2],args[3],args[4],
                                    decoratorString,
                                    whitespaceString
                            }
                    );
                }
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

    @Test(dataProvider = "provideForTestVersionRegex")
    public void testVersionRegex(final Integer major,
                                 final Integer minor,
                                 final Integer year,
                                 final Integer month,
                                 final Integer day,
                                 final String decorator,
                                 final String leadingWhitespace ) {

        // Construct the string:
        final String versionString = String.format(
                "%s%s%d.%d.%4d%02d%02d%s",
                 DataSourceUtils.MANIFEST_VERSION_LINE_START,leadingWhitespace,
                 major,minor,year,month,day,decorator
        );

        final Matcher matcher = DataSourceUtils.VERSION_PATTERN.matcher(versionString);

        Assert.assertTrue(matcher.matches());

        final Integer versionMajor     = Integer.valueOf(matcher.group(1));
        final Integer versionMinor     = Integer.valueOf(matcher.group(2));
        final Integer versionYear      = Integer.valueOf(matcher.group(3));
        final Integer versionMonth     = Integer.valueOf(matcher.group(4));
        final Integer versionDay       = Integer.valueOf(matcher.group(5));
        final String  versionDecorator = matcher.group(6);

        Assert.assertEquals( versionMajor, major );
        Assert.assertEquals( versionMinor, minor );
        Assert.assertEquals( versionYear, year );
        Assert.assertEquals( versionMonth, month );
        Assert.assertEquals( versionDay, day );
        Assert.assertEquals( versionDecorator, decorator );
    }

}
