package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import org.aeonbits.owner.util.Collections;
import org.apache.commons.lang.RandomStringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.GregorianCalendar;
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
        baseArgs.add(DataSourceUtils.MIN_DATE.get(Calendar.YEAR));
        baseArgs.add(DataSourceUtils.MIN_DATE.get(Calendar.MONTH));
        baseArgs.add(DataSourceUtils.MIN_DATE.get(Calendar.DAY_OF_MONTH));

        final ArrayList<Integer> goodRange = new ArrayList<>();
        goodRange.add(DataSourceUtils.MAX_MAJOR_VERSION_NUMBER - DataSourceUtils.MIN_MAJOR_VERSION_NUMBER);
        goodRange.add(DataSourceUtils.MAX_MINOR_VERSION_NUMBER - DataSourceUtils.MIN_MINOR_VERSION_NUMBER);
        goodRange.add(
                DataSourceUtils.MAX_DATE.get(Calendar.YEAR) - DataSourceUtils.MIN_DATE.get(Calendar.YEAR)
        );
        goodRange.add(
                DataSourceUtils.MAX_DATE.get(Calendar.MONTH) - DataSourceUtils.MIN_DATE.get(Calendar.MONTH)
        );
        goodRange.add(
                DataSourceUtils.MAX_DATE.get(Calendar.DAY_OF_MONTH) - DataSourceUtils.MIN_DATE.get(Calendar.DAY_OF_MONTH)
        );

        final int MAX_OFFSET = 5;

        for ( int i = 0; i < baseArgs.size(); ++i ) {
            for ( int offset = -MAX_OFFSET; offset < goodRange.get(i) + MAX_OFFSET; ++offset){

                final ArrayList<Object> argList = new ArrayList<>(baseArgs.subList(0, i));
                argList.add(baseArgs.get(i) + offset);

                if ( i < baseArgs.size() - 1 ) {
                    argList.addAll(baseArgs.subList(i + 1, baseArgs.size()));
                }

                try {
                    final Calendar releaseDate = new GregorianCalendar(
                            (int) argList.get(2), Utils.getCalendarMonth((int) argList.get(3)), (int) argList.get(4)
                    );

                    boolean passes = ((Integer) argList.get(0) >= DataSourceUtils.MIN_MAJOR_VERSION_NUMBER) && ((Integer) argList.get(0) <= DataSourceUtils.MAX_MAJOR_VERSION_NUMBER);
                    passes = passes && ((Integer) argList.get(1) >= DataSourceUtils.MIN_MINOR_VERSION_NUMBER) && ((Integer) argList.get(1) <= DataSourceUtils.MAX_MINOR_VERSION_NUMBER);
                    passes = passes &&
                            (!releaseDate.before(DataSourceUtils.MIN_DATE)) &&
                            (!releaseDate.after(DataSourceUtils.MAX_DATE));

                    testArgs.add(
                            Collections.list(
                                    argList.get(0),
                                    argList.get(1),
                                    releaseDate,
                                    passes
                            ).toArray()
                    );
                }
                catch (final IllegalArgumentException ex) {
                    // If this happened we gave the code an invalid month.  We should ignore this case.
                }
            }
        }

        // Tests for the min version:
        final ArrayList<Integer> minArgs = new ArrayList<>();
        minArgs.add(DataSourceUtils.MIN_MAJOR_VERSION_NUMBER);
        minArgs.add(DataSourceUtils.MIN_MINOR_VERSION_NUMBER);
        minArgs.add(DataSourceUtils.MIN_DATE.get(Calendar.YEAR));
        minArgs.add(DataSourceUtils.MIN_DATE.get(Calendar.MONTH));
        minArgs.add(DataSourceUtils.MIN_DATE.get(Calendar.DAY_OF_MONTH));

        for ( int offset = -2 ; offset < -1; ++offset ) {
            for ( int i = 0; i < minArgs.size(); ++i ) {

                final ArrayList<Object> argList = new ArrayList<>(minArgs.subList(0, i));
                argList.add(minArgs.get(i) + offset);

                if ( i < minArgs.size() - 1) {
                    argList.addAll(minArgs.subList(i + 1, minArgs.size()));
                }

                try {
                    testArgs.add(
                            Collections.list(
                                    argList.get(0),
                                    argList.get(1),
                                    new GregorianCalendar(
                                            (int) argList.get(2), (int) argList.get(3), (int) argList.get(4)
                                    ),
                                    false
                            ).toArray()
                    );
                }
                catch (final IllegalArgumentException ex) {
                    // If this happened we gave the code an invalid month.  We should ignore this case.
                }
            }
        }

        // Tests for max version:
        final ArrayList<Integer> maxArgs = new ArrayList<>();
        maxArgs.add(DataSourceUtils.MAX_MAJOR_VERSION_NUMBER);
        maxArgs.add(DataSourceUtils.MAX_MINOR_VERSION_NUMBER);
        maxArgs.add(DataSourceUtils.MAX_DATE.get(Calendar.YEAR));
        maxArgs.add(DataSourceUtils.MAX_DATE.get(Calendar.MONTH));
        maxArgs.add(DataSourceUtils.MAX_DATE.get(Calendar.DAY_OF_MONTH));

        for ( int offset = 1 ; offset < 3; ++offset ) {
            for ( int i = 0; i < maxArgs.size(); ++i ) {

                final ArrayList<Object> argList = new ArrayList<>(maxArgs.subList(0, i));
                argList.add(maxArgs.get(i) + offset);

                if ( i < maxArgs.size() - 1) {
                    argList.addAll(maxArgs.subList(i + 1, maxArgs.size()));
                }

                try {
                    testArgs.add(
                        Collections.list(
                            argList.get(0),
                            argList.get(1),
                            new GregorianCalendar(
                                    (int)argList.get(2), (int) argList.get(3), (int)argList.get(4)
                            ),
                            false
                        ).toArray()
                    );
                }
                catch (final IllegalArgumentException ex) {
                    // If this happened we gave the code an invalid month.  We should ignore this case.
                }
            }
        }

        // Some specific test cases to prevent regression of version checks:
        // 1 Month after OK release date, but 1 day before OK release date (should pass):
        final Calendar c1 = new GregorianCalendar(DataSourceUtils.MIN_DATE.get(Calendar.YEAR), DataSourceUtils.MIN_DATE.get(Calendar.MONTH) + 1, DataSourceUtils.MIN_DATE.get(Calendar.DAY_OF_MONTH) - 1);
        testArgs.add( new Object[] { DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, c1, true } );
        // 1 Year after OK release date, but 1 day before OK release date (should pass):
        final Calendar c2 = new GregorianCalendar(DataSourceUtils.MIN_DATE.get(Calendar.YEAR)+1, DataSourceUtils.MIN_DATE.get(Calendar.MONTH) , DataSourceUtils.MIN_DATE.get(Calendar.DAY_OF_MONTH) - 1);
        testArgs.add( new Object[] { DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, c2, true } );
        // 1 Year after OK release date, but 1 month before OK release date (should pass):
        final Calendar c3 = new GregorianCalendar(DataSourceUtils.MIN_DATE.get(Calendar.YEAR)+1, DataSourceUtils.MIN_DATE.get(Calendar.MONTH) - 1, DataSourceUtils.MIN_DATE.get(Calendar.DAY_OF_MONTH));
        testArgs.add( new Object[] { DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, c3, true } );
        // 1 Year after OK release date, but 1 month and 1 day before OK release date (should pass):
        final Calendar c4 = new GregorianCalendar(DataSourceUtils.MIN_DATE.get(Calendar.YEAR)+1, DataSourceUtils.MIN_DATE.get(Calendar.MONTH) - 1, DataSourceUtils.MIN_DATE.get(Calendar.DAY_OF_MONTH) - 1);
        testArgs.add( new Object[] { DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, c4, true } );

        // A value in the middle:
        final Calendar c5 = new GregorianCalendar(DataSourceUtils.MIN_DATE.get(Calendar.YEAR)+1, DataSourceUtils.MIN_DATE.get(Calendar.MONTH), DataSourceUtils.MIN_DATE.get(Calendar.DAY_OF_MONTH));
        testArgs.add( new Object[] { DataSourceUtils.MIN_MAJOR_VERSION_NUMBER, DataSourceUtils.MIN_MINOR_VERSION_NUMBER, c5, true } );

        // 1 Month before OK max date, but 1 day before OK max release date (should pass):
        final Calendar c6 = new GregorianCalendar(DataSourceUtils.MAX_DATE.get(Calendar.YEAR), DataSourceUtils.MAX_DATE.get(Calendar.MONTH) - 1, DataSourceUtils.MAX_DATE.get(Calendar.DAY_OF_MONTH) + 1);
        testArgs.add( new Object[] { DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, DataSourceUtils.MAX_MINOR_VERSION_NUMBER, c6, true } );
        // 1 Year before OK max release date, but 1 day before OK max release date (should pass):
        final Calendar c7 = new GregorianCalendar(DataSourceUtils.MAX_DATE.get(Calendar.YEAR) - 1  , DataSourceUtils.MAX_DATE.get(Calendar.MONTH), DataSourceUtils.MAX_DATE.get(Calendar.DAY_OF_MONTH) + 1);
        testArgs.add( new Object[] { DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, DataSourceUtils.MAX_MINOR_VERSION_NUMBER, c7, true } );
        // 1 Year before OK max release date, but 1 month before OK max release date (should pass):
        final Calendar c8 = new GregorianCalendar(DataSourceUtils.MAX_DATE.get(Calendar.YEAR) - 1, DataSourceUtils.MAX_DATE.get(Calendar.MONTH) - 1, DataSourceUtils.MAX_DATE.get(Calendar.DAY_OF_MONTH));
        testArgs.add( new Object[] { DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, DataSourceUtils.MAX_MINOR_VERSION_NUMBER, c8, true } );
        // 1 Year before OK max release date, but 1 month and 1 day before OK max release date (should pass):
        final Calendar c9 = new GregorianCalendar(DataSourceUtils.MAX_DATE.get(Calendar.YEAR)-1, DataSourceUtils.MAX_DATE.get(Calendar.MONTH) - 1, DataSourceUtils.MAX_DATE.get(Calendar.DAY_OF_MONTH) + 1);
        testArgs.add( new Object[] { DataSourceUtils.MAX_MAJOR_VERSION_NUMBER, DataSourceUtils.MAX_MINOR_VERSION_NUMBER, c9, true } );

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

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForValidateVersionInformation")
    public void testValidateVersionInformation(final Integer major,
                                               final Integer minor,
                                               final Calendar releaseDate,
                                               final Boolean expected) {
        Assert.assertEquals(
                DataSourceUtils.validateVersionInformation(
                    major,
                    minor,
                    releaseDate.get(Calendar.YEAR),
                    releaseDate.get(Calendar.MONTH)+1,
                    releaseDate.get(Calendar.DAY_OF_MONTH)),
            expected.booleanValue()
        );
    }

    @Test(dataProvider = "provideForTestVersionRegex")
    public void testVersionRegex(final Integer major,
                                 final Integer minor,
                                 final Calendar releaseDate,
                                 final String decorator,
                                 final String leadingWhitespace ) {

        // Construct the string:
        final String versionString = String.format(
                "%s%s%d.%d.%4d%02d%02d%s",
                DataSourceUtils.MANIFEST_VERSION_LINE_START,
                leadingWhitespace,
                major,
                minor,
                releaseDate.get(Calendar.YEAR),
                releaseDate.get(Calendar.MONTH) + 1,
                releaseDate.get(Calendar.DAY_OF_MONTH),
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
        Assert.assertEquals( versionYear.intValue(), releaseDate.get(Calendar.YEAR) );
        Assert.assertEquals( versionMonth.intValue(), releaseDate.get(Calendar.MONTH) + 1 );
        Assert.assertEquals( versionDay.intValue(), releaseDate.get(Calendar.DAY_OF_MONTH) );
        Assert.assertEquals( versionDecorator, decorator );
    }

}
