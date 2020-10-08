package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import org.apache.commons.lang.RandomStringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.time.LocalDate;
import java.time.Month;
import java.util.ArrayList;
import java.util.Arrays;
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

        // =============================================================
        // First do tests that check the major / minor version numbers:
        // -------------------------------------------------------------

        final int MINOR_IT_MAX_VAL = 20;

        // major < MIN_MAJOR (minor shouldn't matter)
        for (int minor = 0; minor < MINOR_IT_MAX_VAL; ++minor) {
            testArgs.add(
                    new Object[]{
                            DataSourceUtils.MIN_MAJOR_VERSION_NUMBER - 1, minor, DataSourceUtils.MIN_DATE, false
                    }
            );
        }

        // major == MIN_MAJOR (minor does matter)
        for (int minor = 0; minor < MINOR_IT_MAX_VAL; ++minor) {
            testArgs.add(
                    new Object[]{
                            DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, minor, DataSourceUtils.MIN_DATE,
                            (minor == DataSourceUtils.MIN_MINOR_VERSION_NUMBER || minor == DataSourceUtils.MAX_MINOR_VERSION_NUMBER)
                    }
            );
        }

        // major > MIN_MAJOR and < MAX_MAJOR (minor shouldn't matter)
        // Can't currently test here - version numbers use integer values and right now
        // MIN_MAJOR_VERSION_NUMBER == MAX_MAJOR_VERSION_NUMBER
        // TODO: Add in test when possible.

        // major == MAX_MAJOR (minor does matter)
        for (int minor = 0; minor < MINOR_IT_MAX_VAL; ++minor) {
            testArgs.add(
                    new Object[]{
                            DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, minor, DataSourceUtils.MAX_DATE,
                            (minor == DataSourceUtils.MIN_MINOR_VERSION_NUMBER || minor == DataSourceUtils.MAX_MINOR_VERSION_NUMBER)
                    }
            );
        }

        // major > MAX_MAJOR (minor shouldn't matter)
        for (int minor = 0; minor < MINOR_IT_MAX_VAL; ++minor) {
            testArgs.add(
                    new Object[]{
                            DataSourceUtils.MAX_MAJOR_VERSION_NUMBER + 1, minor, DataSourceUtils.MAX_DATE, false
                    }
            );
        }

        // =============================================================
        // Next do tests that check the dates themselves:
        // -------------------------------------------------------------

        // ----------------------------------------
        // Dates in range:
        // ----------------------------------------

        // Min version with minimally / maximally acceptable dates:
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_DATE,
                        true
                }
        );
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_DATE,
                        true
                }
        );

        // Max version with minimally / maximally acceptable dates:
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_DATE,
                        true
                }
        );
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_DATE,
                        true
                }
        );

        // Year / Month inside acceptable dates with any valid day:
        for (int day = 1; day < 30; ++day) {
            testArgs.add(
                    new Object[]{
                            DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                            DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                            LocalDate.of(2020, 4, day),
                            true
                    }
            );
        }

        // 1 day inside the acceptable release window:
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MIN_DATE.getYear(),
                                DataSourceUtils.MIN_DATE.getMonth(),
                                DataSourceUtils.MIN_DATE.getDayOfMonth() + 1
                        ),
                        true
                }
        );
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MAX_DATE.getYear(),
                                DataSourceUtils.MAX_DATE.getMonth(),
                                DataSourceUtils.MAX_DATE.getDayOfMonth() - 1
                        ),
                        true
                }
        );

        // ----------------------------------------
        // Dates out of range:
        // ----------------------------------------

        // 1 day outside the acceptable release window:
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MIN_DATE.getYear(),
                                DataSourceUtils.MIN_DATE.getMonth(),
                                DataSourceUtils.MIN_DATE.getDayOfMonth() - 1
                        ),
                        false
                }
        );
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MAX_DATE.getYear(),
                                DataSourceUtils.MAX_DATE.getMonth(),
                                DataSourceUtils.MAX_DATE.getDayOfMonth() + 1
                        ),
                        false
                }
        );

        // 1 month outside the acceptable release window:

        // Note: Min date month is January, so cannot test one month before it -
        //       that test case is covered by the year being below the minimum.
        // TODO: Add in test when possible.

        testArgs.add(
                new Object[]{
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MAX_DATE.getYear(),
                                DataSourceUtils.MAX_DATE.getMonthValue() + 1,
                                DataSourceUtils.MAX_DATE.getDayOfMonth()
                        ),
                        false
                }
        );

        // 1 year outside the acceptable release window:
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MIN_DATE.getYear() - 1,
                                DataSourceUtils.MIN_DATE.getMonth(),
                                DataSourceUtils.MIN_DATE.getDayOfMonth()
                        ),
                        false
                }
        );
        testArgs.add(
                new Object[]{
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        LocalDate.of(
                                DataSourceUtils.MAX_DATE.getYear() + 1,
                                DataSourceUtils.MAX_DATE.getMonth(),
                                DataSourceUtils.MAX_DATE.getDayOfMonth()
                        ),
                        false
                }
        );

        // Valid Year, month outside bounds, all day values:
        for (int day = 1; day < 30; ++day) {
            testArgs.add(
                    new Object[]{
                            DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                            DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                            LocalDate.of(2020, 9, day),
                            false
                    }
            );
        }

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

                    // Some sanity checks here for proper version numbers:
                    if (((Integer)args[0]) < 0 || ((Integer)args[1]) < 0){
                        continue;
                    }

                    final String whitespaceString = whitespace != 0 ? "\t \t \t " : " ";
                    final String decoratorString  = decoratorCount  !=0 ? RandomStringUtils.randomAlphanumeric(decoratorCount) : "";

                    testArgs.add(
                            new Object[] {
                                    args[0],args[1],args[2],
                                    decoratorString,
                                    whitespaceString
                            }
                    );
                }
            }
        }

        return testArgs.iterator();

    }

    @DataProvider
    private Object[][] provideForGetDataSourceVersionString() {
        return new Object[][] {
                {
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_DATE,
                        "v1.6.20190124"
                },
                {
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_DATE,
                        "v1.7.20200521"
                },
                {
                        52,
                        13,
                        LocalDate.of(71999, Month.DECEMBER, 7),
                        "v52.13.719991207"
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForGetDataSourceVersionString")
    public void testGetDataSourceVersionString(final int major,
                                               final int minor,
                                               final LocalDate date,
                                               final String expected) {
        Assert.assertEquals(
                DataSourceUtils.getDataSourceVersionString(major, minor, date),
                expected
        );
    }

    @Test(dataProvider = "provideForValidateVersionInformation")
    public void testValidateVersionInformation(final Integer major,
                                               final Integer minor,
                                               final LocalDate releaseDate,
                                               final Boolean expected) {
        Assert.assertEquals(
                DataSourceUtils.validateVersionInformation(
                    major,
                    minor,
                    releaseDate.getYear(),
                    releaseDate.getMonthValue(),
                    releaseDate.getDayOfMonth()
                ),
            expected.booleanValue()
        );
    }

    @Test
    public void testGetDataSourceMaxVersionString() {
        Assert.assertEquals(
                DataSourceUtils.getDataSourceMaxVersionString(),
                DataSourceUtils.getDataSourceVersionString(
                        DataSourceUtils.MAX_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MAX_DATE
                )
        );
    }

    @Test
    public void testGetDataSourceMinVersionString() {
        Assert.assertEquals(
                DataSourceUtils.getDataSourceMinVersionString(),
                DataSourceUtils.getDataSourceVersionString(
                        DataSourceUtils.MIN_MAJOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_MINOR_VERSION_NUMBER,
                        DataSourceUtils.MIN_DATE
                )
        );
    }

    @Test(dataProvider = "provideForTestVersionRegex")
    public void testVersionRegex(final Integer major,
                                 final Integer minor,
                                 final LocalDate releaseDate,
                                 final String decorator,
                                 final String leadingWhitespace ) {

        // Construct the string:
        final String versionString = String.format(
                "%s%s%d.%d.%4d%02d%02d%s",
                DataSourceUtils.MANIFEST_VERSION_LINE_START,
                leadingWhitespace,
                major,
                minor,
                releaseDate.getYear(),
                releaseDate.getMonthValue(),
                releaseDate.getDayOfMonth(),
                decorator
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
        Assert.assertEquals( versionYear.intValue(), releaseDate.getYear() );
        Assert.assertEquals( versionMonth.intValue(), releaseDate.getMonthValue() );
        Assert.assertEquals( versionDay.intValue(), releaseDate.getDayOfMonth() );
        Assert.assertEquals( versionDecorator, decorator );
    }

    @Test(groups = {"bucket", "cloud"})
    public void testCurrentDataSourcesAvailable() {

        // Make sure that we can actually access the current min and max versions of the data sources.
        final List<String> dataSourcesBasePaths = Arrays.asList(
                DataSourceUtils.DATA_SOURCES_BUCKET_PATH +
                        DataSourceUtils.DATA_SOURCES_NAME_PREFIX +
                        "." + DataSourceUtils.getDataSourceMinVersionString(),
                DataSourceUtils.DATA_SOURCES_BUCKET_PATH +
                        DataSourceUtils.DATA_SOURCES_NAME_PREFIX +
                        "." + DataSourceUtils.getDataSourceMaxVersionString()
        );
        for (final String basePath : dataSourcesBasePaths) {
            for (final String useCaseModifier : Arrays.asList(DataSourceUtils.DS_SOMATIC_NAME_MODIFIER, DataSourceUtils.DS_GERMLINE_NAME_MODIFIER)) {
                // Data Sources File:
                IOUtils.assertFileIsReadable(
                        IOUtils.getPath(
                            basePath + useCaseModifier + DataSourceUtils.DS_EXTENSION
                        )
                );

                // Checksum for file:
                IOUtils.assertFileIsReadable(
                        IOUtils.getPath(
                                basePath + useCaseModifier + DataSourceUtils.DS_CHECKSUM_EXTENSION
                        )
                );
            }
        }
    }

}
