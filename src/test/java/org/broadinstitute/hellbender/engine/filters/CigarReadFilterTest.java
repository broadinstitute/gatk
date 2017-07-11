package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.ReadClipperTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
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
                {"^", false},                           // Start anchor with no content
                {"$", false},                           // End anchor with no content
                {"SHS", false},                         // Soft clips before hard clips
                {"SHSH", false},                        // Soft clips before hard clips
                {"329", false},                         // Number specified without content
                {"M*", false},                          // * is not a wildcard - it is reserved for "unavailable"
                {"*M", false},                          // * is not a wildcard - it is reserved for "unavailable"
                {"<=P", false},                         // Floating less than
                {"10P<>34X", false},                    // Floating less than
                {"10P>>=34X", false},                   // Floating greater than
                {"12H<3SH10P>=34X3S2H", false},         // Numerically specified hard clips inside cigar
        };
    }

    @DataProvider
    public Object[][] createValidCigarFilterStrings() {
        return new Object[][]{
                {"*", true},                            // Unavailable cigar string
                {"H", true},                            // Single hard clip
                {"HS", true},                           // Hard clip, then a soft clip
                {"M%", true},                           // An Alignment Match then a wildcard
                {"HHHHSSSSDSSSSHHH", true},             // Valid clipping bookends and a Deletion
                {"SS==3I<=24DMMMDMDMDDSHHHH", true},    // Valid clipping bookends with Sequence Matches, Insertions, Numerically specified Deletions and Alignment Matches / Deletions
                {"SH", true},                           // Soft clipping then hard clipping
                {"=H", true},                           // Sequence match then hard clipping
                {"105H", true},                         // Numerically specified hard clipping
                {"329H3S", true},                       // Numerically specified Hard clipping then Numerically specified soft clipping
                {"^M$", true},                          // Anchors and a valid CigarElement
                {"M", true},                            // A valid CigarElement
                {"10P", true},                          // A numerically qualified valid CigarElement
                {"=10P", true},                         // A sequence match, then a numerically qualified valid CigarElement
                {"<=9=", true},                         // A numerical-comparator qualified Sequence Match
                {"10P=34X", true},                      // Numerical qualifiers, sequence match
                {"10P<34X", true},                      // Numerical comparator/qualifiers with valid CigarElements
                {"10P<=34X", true},                     // Numerical comparator/qualifiers with valid CigarElements
                {"10P>34X", true},                      // Numerical comparator/qualifiers with valid CigarElements
                {"10P>=34X", true},                     // Numerical comparator/qualifiers with valid CigarElements
                {"10P>=34X3S2H", true},                 // Numerical comparator/qualifiers with valid CigarElements and clipping
                {"12H<3S10P>=34X3S2H", true},           // Numerical comparator/qualifiers with valid CigarElements and clipping
                {"^12H<3S10P>=34X3S2H$", true},         // Numerical comparator/qualifiers with valid CigarElements, clipping, and anchors
        };
    }

    @DataProvider
    public Object[][] createExhaustiveValidCigarList() {

        ArrayList<String> masterList = new ArrayList<>();

        // NOTE: This can be adjusted for larger scale tests,
        // but it takes a LONG time (as you might expect).
        final int MIN_TEST_LENGTH = 1;
        final int MAX_TEST_LENGTH = 5;

        // Generate the tests:
        for (int i = MIN_TEST_LENGTH; i <= MAX_TEST_LENGTH; ++i) {
            ReadClipperTestUtils.generateCigarList(i).forEach(c -> masterList.add(c.toString()));
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

//        "(\\^?(?:[<>]|[<>]=)?\\d*[SHMIDNPX=%]\\$?)"
        return new Object[][] {
                {"^", false},
                {"$", false},
                {"M$", true},
                {"^*$", false},
                {"^%$", true},
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

//        "(\\^?(?:[<>]\\d+|[<>]=\\d+|\\d+)?[SHMIDNPX=%]\\$?)"

        // Create a list of basic match elements for the core of each match element:
        ArrayList<String> basicMatchElements =
                new ArrayList<>(Arrays.asList("S","H","M","I","D","N","P","X","=",
                        Character.toString( CigarReadFilter.CigarMatchElement.WILDCARD )));

        // Create a list of decorators for each match element:
        ArrayList<String> decorators =
                new ArrayList<>(Arrays.asList("", "<", ">", "<=",">="));

        // Create a place to put the combined elements and decorators:
        Object[][] decoratedMatchElements = new Object[ basicMatchElements.size() * (4 * (decorators.size()+1)) ][1];
        int dIndex = 0;

        // Create and insert the new decorated elements into the output array:
        for ( String e : basicMatchElements ) {

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

    @DataProvider
    public Object[][] createSetsOfSingleMatchStrings() {

        // Create a list of basic match elements for the core of each match element:
        ArrayList<String> basicMatchElements =
                new ArrayList<>(Arrays.asList("S","H","M","I","D","N","P","X","=",
                        Character.toString( CigarReadFilter.CigarMatchElement.WILDCARD )));

        // Create a list of decorators for each match element:
        ArrayList<String> decorators =
                new ArrayList<>(Arrays.asList("", "<", ">", "<=",">="));

        // Create an output list:
        ArrayList<ArrayList<String>> outList = new ArrayList<>();

        // Add in the unavailable string to the outlist:
        ArrayList<String> unavailableString = new ArrayList<>();
        unavailableString.add("*");
        outList.add(unavailableString);

        // NOTE: These two settings are specifically set to balance the need for tests but
        //       also speed of tests.
        // The minimum number of elements to put as a single test:
        final int MIN_NUM_ELEMENT_STRINGS = 1;

        // The maximum number of elements to put as a single test:
        final int MAX_NUM_ELEMENT_STRINGS = 4;

        for ( int elLength = MIN_NUM_ELEMENT_STRINGS; elLength < MAX_NUM_ELEMENT_STRINGS; ++elLength ) {

            // Create the permutations:
            Collection<ArrayList<String>> elementPermutations = getAllPermutations(elLength, basicMatchElements);

            // For each one, add decorators as appropriate:
            for (ArrayList<String> permutation : elementPermutations ) {

                // Add trivial cases first:
                outList.add(permutation);

                ArrayList<String> permutationStartAnchor = new ArrayList<>(permutation);
                permutationStartAnchor.set(0, "^" + permutationStartAnchor.get(0));
                outList.add(permutationStartAnchor);

                ArrayList<String> permutationEndAnchor = new ArrayList<>(permutation);
                permutationEndAnchor.set(permutationEndAnchor.size()-1, permutationEndAnchor.get(permutationEndAnchor.size()-1) + "$");
                outList.add(permutationEndAnchor);

                ArrayList<String> permutationStartEndAnchor = new ArrayList<>(permutation);
                permutationStartEndAnchor.set(0, "^" + permutationStartEndAnchor.get(0));
                permutationStartEndAnchor.set(permutationStartEndAnchor.size()-1, permutationStartEndAnchor.get(permutationStartEndAnchor.size()-1) + "$");
                outList.add(permutationStartEndAnchor);

                // Now go through the decorators and add them to each element
                // We do this to make it a little faster than doing the combined
                // permutations of the Cigar operators and the decorators:
                for( String d : decorators ) {
                    ArrayList<String> decoratedPermut = new ArrayList<>(permutation);

                    // Decorate the elements:
                    for (int i = 0 ; i < decoratedPermut.size() ; ++i) {
                        Integer num = ThreadLocalRandom.current().nextInt(1, 1000 + 1);
                        decoratedPermut.set(i, d + num.toString() + decoratedPermut.get(i));
                    }

                    // Add in the decorated cases:
                    outList.add(decoratedPermut);

                    ArrayList<String> decoratedPermutStartAnchor = new ArrayList<>(decoratedPermut);
                    decoratedPermutStartAnchor.set(0, "^" + decoratedPermutStartAnchor.get(0));
                    outList.add(decoratedPermutStartAnchor);

                    ArrayList<String> decoratedPermutEndAnchor = new ArrayList<>(decoratedPermut);
                    decoratedPermutEndAnchor.set(decoratedPermutEndAnchor.size()-1, decoratedPermutEndAnchor.get(decoratedPermutEndAnchor.size()-1) + "$");
                    outList.add(decoratedPermutEndAnchor);

                    ArrayList<String> decoratedPermutStartEndAnchor = new ArrayList<>(decoratedPermut);
                    decoratedPermutStartEndAnchor.set(0, "^" + decoratedPermutStartEndAnchor.get(0));
                    decoratedPermutStartEndAnchor.set(decoratedPermutStartEndAnchor.size()-1, decoratedPermutStartEndAnchor.get(decoratedPermutStartEndAnchor.size()-1) + "$");
                    outList.add(decoratedPermutStartEndAnchor);
                }
            }
        }

        // Package everything for the output format:
        Object[][] objectsOut = new Object[outList.size()][1];
        for (int i = 0; i < outList.size(); ++i) {
                objectsOut[i] = new Object[] { outList.get(i) };
        }

        return objectsOut;
    }

    @DataProvider
    public Object[][] createCigarStringsAndMatchesForReadFilterTest() {
        //cigarString, cigarMatchString, shouldMatch

//                "1H1S10M10S30H"
//                "1I9H"
//                "1I1S8H"
//                "1S1I1S7"
//
        return new Object[][]{
                {"*", "*", true},                       // Match the unavailable cigar string
                {"*", "M", false},                      // Match the unavailable cigar string
                {"*", "HSMSH", false},                  // Match the unavailable cigar string
                {"*", "10M", false},                    // Match the unavailable cigar string
                {"*", "<=5M", false},                   // Match the unavailable cigar string

                {"1H", "S", false},                     // Pattern doesn't specify the correct CigarElement
                {"1H", "M", false},                     // Pattern doesn't specify the correct CigarElement
                {"1H", "=", false},                     // Pattern doesn't specify the correct CigarElement
                {"1H", "^107I", false},                 // Pattern doesn't specify the correct CigarElement
                {"1H", "H", true},                      // Pattern matches CigarElement exactly
                {"1H", ">=1H", true},                   // Match equals case for inequalities CigarElements
                {"1H", "<=2H", true},                   // Match equals case for inequalities CigarElements
                {"1H", "<2H", true},                    // Pattern matches CigarElement numerically
                {"2H", ">1H", true},                    // Pattern matches CigarElement numerically
                {"2H2H1H", ">4H", true},                // Pattern matches CigarElement numerically
                {"2H", ">2H", false},                   // Too few CigarElements
                {"1H1H", ">2H", false},                 // Too few CigarElements
                {"2H", ">=3H", false},                  // Too few CigarElements
                {"1H1H", ">=3H", false},                // Too few CigarElements
                {"3H", "<3H", false},                   // Too many CigarElements
                {"2H1H", "<3H", false},                 // Too many CigarElements
                {"4H", "<=3H", false},                  // Too many CigarElements
                {"2H2H", "<=3H", false},                // Too many CigarElements

                {"1P1=", "P=M", false},                 // Match longer than cigar

                {"12H3S10P34X3S2H", "P", true},         // Match inside a larger string
                {"12H3S10P34X3S2H", "10P", true},       // Match inside a larger string
                {"12H3S10P34X3S2H", ">9P", true},       // Match inside a larger string
                {"12H3S10P34X3S2H", "<8P", false},      // Match inside a larger string
                {"12H3S10P34X3S2H", "PX", true},        // Match inside a larger string
                {"12H3S10P34X3S2H", ">5PX", true},      // Match inside a larger string
                {"12H3S10P34X3S2H", "P>2X", true},      // Match inside a larger string
                {"12H3S10P34X3S2H", "10P34X", true},    // Match inside a larger string
                {"12H3S10P34X3S2H", ">9P>33X", true},   // Match inside a larger string
                {"12H3S10P34X3S2H", "<8P>2X", false},   // Match inside a larger string
                {"12H3S10P34X3S2H", ">9PX>2S", true},   // Match inside a larger string
                {"12H3S10P34X3S2H", "PS", false},       // Near-Match inside a larger string
                {"12H3S10P34X3S2H", "10P34XH", false},  // Near-Match inside a larger string
                {"12H3S10P34X3S2H", "<8P>2X>5S", false},// Near-Match inside a larger string

                {"1H1S10M10S30H", "^H", true},          // Match starting anchors
                {"1H1S10M10S30H", "^1H", true},         // Match starting anchors
                {"1H1S10M10S30H", "^2H", false},        // Match starting anchors
                {"1H1S10M10S30H", "^HS", true},         // Match starting anchors
                {"1H1S10M10S30H", "^HSM", true},        // Match starting anchors
                {"1H1S10M10S30H", "^HSM", true},        // Match starting anchors
                {"1H1S10M10S30H", "^H2SM", false},      // Match starting anchors
                {"1H1S10M10S30H", "^HS<10M", false},    // Match starting anchors
                {"1H1S10M10S30H", "^HSMSH", true},      // Match starting anchors

                {"1H1S10M10S30H", "H$", true},          // Match ending anchors
                {"1H1S10M10S30H", "1H$", false},        // Match ending anchors - mismatch number
                {"1H1S10M10S30H", "2H$", false},        // Match ending anchors - mismatch number
                {"1H1S10M10S30H", "30H$", true},        // Match ending anchors
                {"1H1S10M10S30H", ">=2H$", true},       // Match ending anchors
                {"1H1S10M10S30H", ">35H$", false},      // Match ending anchors
                {"1H1S10M10S30H", "SH$", true},         // Match ending anchors
                {"1H1S10M10S30H", "HSM$", false},       // Match ending anchors
                {"1H1S10M10S30H", "MSH$", true},        // Match ending anchors
                {"1H1S10M10S30H", "M13SH", false},      // Match ending anchors
                {"1H1S10M10S30H", "M<10MSH", false},    // Match ending anchors
                {"1H1S10M10S30H", "HSMSH$", true},      // Match ending anchors

                {"1H1S10M10S30H", "^HSMSH$", true},     // Match starting and ending anchors
                {"1H1S10M10S30H", "^HSPSH$", false},    // Match starting and ending anchors
                {"1H1S10M10S30H", "^HS10MSH$", true},   // Match starting and ending anchors
                {"1H1S10M10S30H", "^HS10M10S30H$", true},// Match starting and ending anchors
                {"1H1S10M10S30H", "^HSMS<25H$", false}, // Match starting and ending anchors

                {"1H1S1M1I", "<5%", true},              // Any cigar string of length <5 OK
                {"3H", "<5%", true},                    // Any cigar string of length <5 OK
                {"6S", "<5%", false},                   // Any cigar string of length <5 Fail
                {"1H1S1M1I2S", "<5%", false},           // Any cigar string of length <5 Fail
                {"1H1S1M1I1S", "<5%", false},           // Any cigar string of length <5 Fail

                {"*", "%", false},                      // Wildcard should not match unavailable
                {"1H", "%", true},                      // Wildcard should match any valid cigar string
                {"1H1S", "%", true},                    // Wildcard should match any valid cigar string
                {"1H1S1H", "%", true},                  // Wildcard should match any valid cigar string
                {"12H3S10P34X3S2H", "%", true},         // Wildcard should match any valid cigar string
                {"5H", "5%", true},                     // Any cigar string of length 5 OK
                {"5=", "5%", true},                     // Any cigar string of length 5 OK
                {"2S3M", "5%", true},                   // Any cigar string of length 5 OK
                {"1H1S1M1I1S", "5%", true},             // Any cigar string of length 5 OK
                {"4H", "5%", false},                    // Any cigar string of length 5 Fail
                {"4=", "5%", false},                    // Any cigar string of length 5 Fail
                {"2S1M1I", "5%", false},                // Any cigar string of length 5 Fail
                {"1S2M", "5%", false},                  // Any cigar string of length 5 Fail
                {"1H1S1M1I", "5%", false},              // Any cigar string of length 5 Fail
                {"1M", "M%", false},                    // Alignment match followed by anything
                {"1M1I", "M%", true},                   // Alignment match followed by anything
                {"2M", "M%", true},                     // Alignment match followed by anything
                {"1M1S", "M%", true},                   // Alignment match followed by anything
                {"1S1M", "M%", false},                  // Alignment match followed by anything
                {"1M", "%M", false},                    // Anything followed by alignment match
                {"1I1M", "%M", true},                   // Anything followed by alignment match
                {"2M", "%M", true},                     // Anything followed by alignment match
                {"1S1M", "%M", true},                   // Anything followed by alignment match
                {"1S1M", "%M", true},                   // Anything followed by alignment match
                {"1S1M", "%%", true},                   // Any two cigar elements
                {"1S", "%%", false},                    // Any two cigar elements
                {"1S1M1H", "%%", true},                 // Any two cigar elements


//                {"HS", true},
//                {"M%", true},
//                {"HHHHSSSSDSSSSHHH", true},
//                {"SS==3I<=24DMMMDMDMDDSHHHH", true},
//                {"SH", true},
//                {"=H", true},
//                {"105H", true},
//                {"329H3S", true},
//                {"^M$", true},
//                {"M", true},
//                {"10P", true},
//                {"=10P", true},
//                {"<=9=", true},
//                {"10P=34X", true},
//                {"10P<34X", true},
//                {"10P<=34X", true},
//                {"10P>34X", true},
//                {"10P>=34X", true},
//                {"10P>=34X3S2H", true},
//                {"12H<3S10P>=34X3S2H", true},
//                {"^12H<3S10P>=34X3S2H$", true},
        };
    }

    //==================================================================================================================

    private boolean cigarMatchElementEqualsString(final CigarReadFilter.CigarMatchElement e, final String matchString) {

        boolean result = true;

        // Check the match element against what we expect output to be:

        // Ensure that we have at least some content in the string:
        if ( matchString.length() == 0 ) {
            return false;
        }
        // Check against the unavailable special case:
        else if ( matchString.compareTo("*") == 0) {
            result = result && e.isUnavailable();

            result = result && !e.isAnchoredStart();
            result = result && !e.isAnchoredEnd();
            result = result && !e.isLessThan();
            result = result && !e.isLessThanEqualTo();
            result = result && !e.isGreaterThan();
            result = result && !e.isGreaterThanEqualTo();
            result = result && !e.isWildCard();
            result = result && (e.getLength() == -1);
            result = result && (e.getOperator() == null);
        }
        else {

            // Get the basic operator on which to test:
            char operator;
            operator = matchString.charAt(matchString.length() - 1);

            // Check anchoring:
            if (matchString.contains("^")) {
                result = result && e.isAnchoredStart();
            }
            if (matchString.contains("$")) {
                result = result && e.isAnchoredEnd();

                // We must get the operator again, because we greedily got the "$" instead:
                operator = matchString.charAt(matchString.length() - 2);
            }

            // Check numerical qualifiers:
            if (matchString.contains("<=")) {
                result = result && e.isLessThanEqualTo();
                result = result && !e.isLessThan();
                result = result && !e.isGreaterThan();
                result = result && !e.isGreaterThanEqualTo();
            } else if (matchString.contains("<")) {
                result = result && e.isLessThan();
                result = result && !e.isLessThanEqualTo();
                result = result && !e.isGreaterThan();
                result = result && !e.isGreaterThanEqualTo();
            } else if (matchString.contains(">=")) {
                result = result && e.isGreaterThanEqualTo();
                result = result && !e.isLessThan();
                result = result && !e.isLessThanEqualTo();
                result = result && !e.isGreaterThan();
            } else if (matchString.contains(">")) {
                result = result && e.isGreaterThan();
                result = result && !e.isLessThan();
                result = result && !e.isLessThanEqualTo();
                result = result && !e.isGreaterThanEqualTo();
            }

            // Make sure the operators are the same:
            if (matchString.contains(Character.toString( CigarReadFilter.CigarMatchElement.WILDCARD ))) {
                result = result && e.isWildCard();
            } else {
                result = result && (e.getOperator() == CigarOperator.characterToEnum(operator));
            }

            // Now check the numbers:
            Pattern numberPattern = Pattern.compile("(\\d+)");
            if (numberPattern.matcher(matchString).matches()) {

                // The match string contains numbers, we must check the numbers:
                int numberInMatchPattern = Integer.valueOf(numberPattern.matcher(matchString).group(0));

                result = result && (numberInMatchPattern == e.getLength());
            }
        }

        return result;
    }

    /**
     * Creates a {@link Collection} of all permutations of the given {@link Collection} of the given length.
     * @param length Number of elements to include in the permutations.
     * @param objectList Collection from which to create permutations.
     * @return A {@link Collection} of permutations of the given length of c.
     */
    private <T extends Object>
    Collection<ArrayList<T>> getAllPermutations(int length, List<T> objectList) {
        LinkedList<ArrayList<T>> outList = new LinkedList<>();

        // If we want permutations of length zero, there is only one - the empty set:
        if ( length == 0) {
            return outList;
        }

        byte[] elementCombination = new byte[length];

        Arrays.fill(elementCombination, (byte) 0);
        int currentIndex = 0;
        while (true) {

            ArrayList<T> permutation = new ArrayList<>();

            for (byte i : elementCombination) {
                permutation.add(objectList.get(i));
            }

            outList.add( permutation );    // create the cigar

            boolean currentIndexChanged = false;
            while (currentIndex < length && elementCombination[currentIndex] == objectList.size() - 1) {
                currentIndex++;                                            // find the next index to increment
                currentIndexChanged = true;                                // keep track of the fact that we have changed indices!
            }

            if (currentIndex == length)                             // if we hit the end of the array, we're done.
                break;

            elementCombination[currentIndex]++;                              // otherwise advance the current index

            if (currentIndexChanged) {                                     // if we have changed index, then...
                for (int i = 0; i < currentIndex; i++) {
                    elementCombination[i] = 0;
                }                                                           // reset everything from 0->currentIndex
                currentIndex = 0;                                          // go back to the first index
            }
        }

        return outList;
    }

    //==================================================================================================================

    @Test(dataProvider="createCigarFilterStrings")
    public void testValidate(String cigarFilterString, boolean isValid) {
        Assert.assertEquals( CigarReadFilter.validate(cigarFilterString), isValid );
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

        CigarReadFilter cigarReadFilter = new CigarReadFilter();
        CigarReadFilter.CigarMatchElement e = cigarReadFilter.createMatchElementFromSingleMatchString(matchString);

        // Check the match element against what we expect output to be:

        Assert.assertTrue(cigarMatchElementEqualsString(e, matchString));
    }

    //==================================================================================================================

    @Test(dataProvider="createSetsOfSingleMatchStrings")
    public void testExtractMatchElementsFromString( ArrayList<String> matchStrings ) {

        // First concatenate the strings in the array list:
        StringBuilder matchStringBuilder = new StringBuilder();
        for ( String s : matchStrings ) {
            matchStringBuilder.append(s);
        }

        CigarReadFilter cigarReadFilter = new CigarReadFilter();
        List<CigarReadFilter.CigarMatchElement> c =
                cigarReadFilter.extractMatchElementsFromString(matchStringBuilder.toString());

        for ( int i = 0 ; i < matchStrings.size(); ++i ) {
            Assert.assertTrue(cigarMatchElementEqualsString(c.get(i), matchStrings.get(i)));
        }
    }

    @Test(dataProvider="createExhaustiveValidCigarList")
    public void testValidateExhaustivly(String cigarFilterString) {
        Assert.assertEquals( CigarReadFilter.validate(cigarFilterString), true );
    }

    //==================================================================================================================


    @Test(dataProvider="createCigarStringsAndMatchesForReadFilterTest")
    public void testReadFilterTestMethod(String cigarString, String cigarMatchString, boolean shouldMatch) {
        //FIXME

        GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigarString);
        CigarReadFilter cigarReadFilter = new CigarReadFilter(cigarMatchString);

        Assert.assertEquals(cigarReadFilter.test(read), shouldMatch);
    }

}
