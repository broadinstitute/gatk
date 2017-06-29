package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.ReadClipperTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.regex.Pattern;
import java.util.stream.Stream;

public class CigarReadFilterTest {

    @DataProvider
    public Object[][] createGatkReadsForFilterTest() {
        return new Object[][] {
                { ReadClipperTestUtils.makeReadFromCigar("D") }
        };
    }

    @DataProvider
    public Object[][] createInvalidCigarFilterStrings() {
        return new Object[][]{
                {"^", false},
                {"$", false},
                {"SHS", false},
                {"SHSH", false},
                {"329", false},
                {"M%", false},
                {"<=P", false},
                {"10P<>34X", false},
                {"10P<>=34X", false},
                {"12H<3SH10P>=34X3S2H", false},
        };
    }

    @DataProvider
    public Object[][] createValidCigarFilterStrings() {
        return new Object[][]{
                {"*", true},
                {"H", true},
                {"HS", true},
                {"HHHHSSSSDSSSSHHH", true},
                {"SS==3I<=24DMMMDMDMDDSHHHH", true},
                {"SH", true},
                {"=H", true},
                {"105H", true},
                {"329H3S", true},
                {"^M$", true},
                {"M", true},
                {"M*", true},
                {"*M", true},
                {"10P", true},
                {"=10P", true},
                {"<=9=", true},
                {"10P=34X", true},
                {"10P<34X", true},
                {"10P<=34X", true},
                {"10P>34X", true},
                {"10P>=34X", true},
                {"10P>=34X3S2H", true},
                {"12H<3S10P>=34X3S2H", true},
                {"^12H<3S10P>=34X3S2H$", true},
        };
    }

    @DataProvider
    public Object[][] createExhaustiveValidCigarList() {

        ArrayList<String> masterList = new ArrayList<>();

        // NOTE: This can be adjusted for larger scale tests,
        // but it takes a LONG time (as you might expect).
        final int MAX_TEST_LENGTH = 5;

        // Generate the tests:
        for (int i = 1; i <= MAX_TEST_LENGTH; ++i) {
            ReadClipperTestUtils.generateCigarList(i).stream().forEach(c -> masterList.add(c.toString()));
        }

        // Assign the tests to an output object:
        Object[][] exhaustiveValues = new Object[masterList.size()][1];
        for (int i = 0; i < masterList.size(); ++i) {
            exhaustiveValues[i][0] = masterList.get(i);
        }

        return exhaustiveValues;
    }

    @DataProvider
    public Object[][] createCigarFilterStrings() {

        return Stream.concat(   Arrays.stream( createValidCigarFilterStrings() ),
                                Arrays.stream( createInvalidCigarFilterStrings() )
                            ).toArray(Object[][]::new);

    }

    @DataProvider
    public Object[][] createCigarMatchElementStrings() {

//        "(\\^?(?:[<>]|[<>]=)?\\d*[SHMIDNPX=*]\\$?)"
        return new Object[][] {
                {"^", false},
                {"$", false},
                {"M$", true},
                {"^*$", true},
                {"107X$", true},
                {"=19X", false},
                {"<=19I", true},
                {">=19M", true},
                {">19D", true},
                {"<19N", true},
        };
    }

    @DataProvider
    public Object[][] createSingleMatchStrings() {

//        "(\\^?(?:[<>]\\d+|[<>]=\\d+|\\d+)?[SHMIDNPX=*]\\$?)"

        // Create a list of basic match elements for the core of each match element:
        ArrayList<String> basicMatchElements =
                new ArrayList<>(Arrays.asList("S","H","M","I","D","N","P","X","=","*"));

        // Create a list of decorators for each match element:
        ArrayList<String> decorators =
                new ArrayList<>(Arrays.asList("", "<", ">", "<=",">="));

        // Create a place to put the combined elements and decorators:
        Object[][] decoratedMatchElements = new Object[ basicMatchElements.size() * (4 * (decorators.size()+1)) ][1];
        int dIndex = 0;

        // Create and insert the new decorated elements into the output array:
        for (int i = 0 ; i < basicMatchElements.size(); ++i) {

            String e = basicMatchElements.get(i);

            // Add in the simple element with anchors:
            decoratedMatchElements[dIndex++] = new Object[] { e };
            decoratedMatchElements[dIndex++] = new Object[] { e + "$" };
            decoratedMatchElements[dIndex++] = new Object[] { "^" + e };
            decoratedMatchElements[dIndex++] = new Object[] { "^" + e + "$" };

            for (String d : decorators) {

                Integer num = ThreadLocalRandom.current().nextInt(1, 1000 + 1);
                String eDecorated = d + num.toString() + e;
                decoratedMatchElements[dIndex++] = new Object[] { eDecorated };
                decoratedMatchElements[dIndex++] = new Object[] { eDecorated + "$" };
                decoratedMatchElements[dIndex++] = new Object[] { "^" + eDecorated };
                decoratedMatchElements[dIndex++] = new Object[] { "^" + eDecorated + "$" };
            }
        }

        return decoratedMatchElements;
    }

    //==================================================================================================================

    @Test(dataProvider="createCigarFilterStrings")
    public void testValidate(String cigarFilterString, boolean isValid) {
        Assert.assertEquals( CigarReadFilter.validate(cigarFilterString), isValid );
    }

    @Test(dataProvider="createExhaustiveValidCigarList")
    public void testValidateExhaustivly(String cigarFilterString) {
        Assert.assertEquals( CigarReadFilter.validate(cigarFilterString), true );
    }

    @Test(dataProvider="createCigarMatchElementStrings")
    public void testNextCigarMatchElementPattern(String cigarFilterString, boolean isValid) {
        Assert.assertEquals( CigarReadFilter.nextCigarMatchElementPattern.matcher(cigarFilterString).matches(), isValid );
    }

    @Test(dataProvider="createValidCigarFilterStrings")
    public void testSetDescriptionValid(String cigarFilterString, boolean isValid) {
        CigarReadFilter cigarReadFilter = new CigarReadFilter();

        // Set the CigarReadFilter to have the given description:
        cigarReadFilter.setDescription(cigarFilterString);

        // Verify that the CigarReadDescription is the same as that which was given:
        Assert.assertEquals(cigarReadFilter.getDescription(), cigarFilterString);
    }

    @Test(dataProvider="createInvalidCigarFilterStrings", expectedExceptions=UserException.BadInput.class)
    public void testSetDescriptionInvalid(String cigarFilterString, boolean isValid) {
        CigarReadFilter cigarReadFilter = new CigarReadFilter();

        // Set the CigarReadFilter to have the given description:
        // NOTE: This should throw an exception each time
        cigarReadFilter.setDescription(cigarFilterString);
    }

    @Test(dataProvider="createValidCigarFilterStrings")
    public void testConstructorValidDescription(String cigarFilterString, boolean isValid) {
        CigarReadFilter cigarReadFilter = new CigarReadFilter(cigarFilterString);

        // Verify that the CigarReadDescription is the same as that which was given:
        Assert.assertEquals(cigarReadFilter.getDescription(), cigarFilterString);
    }

    @Test(dataProvider="createInvalidCigarFilterStrings", expectedExceptions=UserException.BadInput.class)
    public void testConstructorInvalidDescription(String cigarFilterString, boolean isValid) {
        CigarReadFilter cigarReadFilter = new CigarReadFilter(cigarFilterString);
    }


    @Test(dataProvider="createSingleMatchStrings")
    public void testCreateMatchElementFromSingleMatchString(String matchString) {
        //FIXME
        CigarReadFilter cigarReadFilter = new CigarReadFilter();
        CigarReadFilter.CigarMatchElement e = cigarReadFilter.createMatchElementFromSingleMatchString(matchString);

        // Check the match element against what we expect output to be:
        if (matchString.contains("^")) {
            Assert.assertTrue(e.isAnchoredStart());
        }
        if (matchString.contains("$")) {
            Assert.assertTrue(e.isAnchoredEnd());
        }

        if (matchString.contains("<=")) {
            Assert.assertTrue(e.isLessThanEqualTo());
            Assert.assertFalse(e.isLessThan());
            Assert.assertFalse(e.isGreaterThan());
            Assert.assertFalse(e.isGreaterThanEqualTo());
        }
        else if (matchString.contains("<")) {
            Assert.assertTrue(e.isLessThan());
            Assert.assertFalse(e.isLessThanEqualTo());
            Assert.assertFalse(e.isGreaterThan());
            Assert.assertFalse(e.isGreaterThanEqualTo());
        }
        else if (matchString.contains(">=")) {
            Assert.assertTrue(e.isGreaterThanEqualTo());
            Assert.assertFalse(e.isLessThan());
            Assert.assertFalse(e.isLessThanEqualTo());
            Assert.assertFalse(e.isGreaterThan());
        }
        else if (matchString.contains(">")) {
            Assert.assertTrue(e.isGreaterThan());
            Assert.assertFalse(e.isLessThan());
            Assert.assertFalse(e.isLessThanEqualTo());
            Assert.assertFalse(e.isGreaterThanEqualTo());
        }

        // Get the basic match element on which to test:
        char basicMatchElement;
        if (matchString.contains("$"))
        {
            basicMatchElement = matchString.charAt(matchString.length()-2);
        }
        else {
            basicMatchElement = matchString.charAt(matchString.length()-1);
        }

        // Make sure the operators are the same:
        if (matchString.contains("*")) {
            Assert.assertTrue(e.isWildCard());
        }
        else {
            Assert.assertEquals(e.getOperator(), CigarOperator.characterToEnum(basicMatchElement));
        }

        // Now check the numbers:
        Pattern numberPattern = Pattern.compile("(\\d+)");
        if ( numberPattern.matcher(matchString).matches() ) {

            // The match string contains numbers, we must check the numbers:
            int numberInMatchPattern = Integer.valueOf( numberPattern.matcher(matchString).group(0) );

            Assert.assertEquals(numberInMatchPattern, e.getLength());
        }
    }

//    @Test
//    public void testExtractMatchElementsFromString() {
//        //FIXME
//        throw new RuntimeException("MUST BE IMPLEMENTED");
//    }

//    @Test
//    public void testReadFilterTestMethod() {
//        //FIXME
//        throw new RuntimeException("MUST BE IMPLEMENTED");
//    }

}
