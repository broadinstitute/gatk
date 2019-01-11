package org.broadinstitute.hellbender.testutils;

import htsjdk.samtools.ValidationStringency;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public final class IntegrationTestSpec {

    /** Standard Logger.  */
    protected static final Logger logger = LogManager.getLogger(IntegrationTestSpec.class);

    public static final String DEFAULT_TEMP_EXTENSION = ".tmp";
    public static final String DEFAULT_TEMP_PREFIX = "walktest.tmp_param";

    private final String args;
    private final Class<?> expectedException;
    private final int nOutputFiles;
    private final List<String> expectedFileNames;

    //If this field is set to true, bam files will be compared after they get sorted.
    //This is needed as a workaround because Spark tools don't respect a pre-ordered BAMs
    // and so may create BAMs that are sorted differently than the input (though both orders will be valid).
    private boolean compareBamFilesSorted;

    //Stringency for validation of bams.
    private ValidationStringency validationStringency;

    public IntegrationTestSpec(String args, List<String> expectedFileNames) {
        this.args = args;
        this.nOutputFiles = expectedFileNames.size();
        this.expectedException = null;
        this.expectedFileNames = expectedFileNames;
        this.compareBamFilesSorted = false;
        this.validationStringency = ValidationStringency.DEFAULT_STRINGENCY;
    }

    public IntegrationTestSpec(String args, int nOutputFiles, Class<?> expectedException) {
        if (expectedException == null) {
            throw new IllegalArgumentException("expected exception is null");
        }
        this.args = args;
        this.nOutputFiles = nOutputFiles;
        this.expectedException = expectedException;
        this.expectedFileNames = null;
        this.compareBamFilesSorted = false;
        this.validationStringency = ValidationStringency.DEFAULT_STRINGENCY;
    }

    public boolean expectsException() {
        return expectedException != null;
    }

    public final Class<?> getExpectedException() {
        if (!expectsException())
            throw new GATKException("Tried to get exception for walker test that doesn't expect one");
        return expectedException;
    }

    public void setCompareBamFilesSorted(final boolean compareBamFilesSorted) {
        this.compareBamFilesSorted = compareBamFilesSorted;
    }

    public void setValidationStringency(final ValidationStringency validationStringency) {
        this.validationStringency = validationStringency;
    }

    public String getArgs() {
        return args;
    }

    public Collection<String> expectedFileNames() {
        return expectedFileNames;
    }

    public void executeTest(final String name, CommandLineProgramTester test) throws IOException {
        List<File> tmpFiles = new ArrayList<>();
        for (int i = 0; i < nOutputFiles; i++) {
            String ext = DEFAULT_TEMP_EXTENSION;
            File fl = BaseTest.createTempFile(String.format(DEFAULT_TEMP_PREFIX + ".%d", i), ext);
            tmpFiles.add(fl);
        }

        final String preFormattedArgs = getArgs();
        final String formattedArgs = String.format(preFormattedArgs, tmpFiles.toArray());
        System.out.println(StringUtils.repeat('-', 80));

        if (expectsException()) {
            // this branch handles the case were we are testing that a walker will fail as expected
            executeTest(name, test, null, null, tmpFiles, formattedArgs, getExpectedException());
        } else {
            final List<String> expectedFileNames = new ArrayList<>(expectedFileNames());
            if (!expectedFileNames.isEmpty() && preFormattedArgs.equals(formattedArgs)) {
                throw new GATKException("Incorrect test specification - you're expecting " + expectedFileNames.size() + " file(s) the specified arguments do not contain the same number of \"%s\" placeholders");
            }

            executeTest(name, test, null, expectedFileNames, tmpFiles, formattedArgs, null);
        }
    }

    /**
     * execute the test, given the following:
     *
     * @param testName              the name of the test
     * @param testClass             the object that contains the test
     * @param expectedFileNames     the list of expectedFileNames
     * @param tmpFiles              the temp file corresponding to the expectedFileNames list
     * @param args                  the argument list
     * @param expectedException     the expected exception or null
     */
    private void executeTest(String testName, CommandLineProgramTester testClass, File outputFileLocation, List<String> expectedFileNames, List<File> tmpFiles, String args, Class<?> expectedException) throws IOException {
        if (outputFileLocation != null) {
            args += " -O " + outputFileLocation.getAbsolutePath();
        }
        executeTest(testName, testClass, args, expectedException);

        if (expectedException == null && !expectedFileNames.isEmpty()) {
            assertMatchingFiles(tmpFiles, expectedFileNames, compareBamFilesSorted, validationStringency);
        }
    }

    /**
     * execute the test, given the following:
     *
     * @param testName          the name of the test
     * @param testClass         the object that contains the test
     * @param args              the argument list
     * @param expectedException the expected exception or null
     */
    private void executeTest(String testName, CommandLineProgramTester testClass, String args, Class<?> expectedException) {
        String[] command = Utils.escapeExpressions(args);
        // run the executable
        boolean gotAnException = false;
        try {
            final String now = LocalDateTime.now().format(DateTimeFormatter.ofPattern("HH:mm:ss"));
            System.out.println(String.format("[%s] Executing test %s:%s", now, testClass.getClass().getSimpleName(), testName));
            testClass.runCommandLine(command);
        } catch (Exception e) {
            gotAnException = true;
            if (expectedException == null) {
                // we didn't expect an exception but we got one :-(
                throw new RuntimeException(e);
            }
            // we expect an exception
            if (!expectedException.isInstance(e)) {
                final String message = String.format("Test %s:%s expected exception %s but instead got %s with error message %s",
                        testClass, testName, expectedException, e.getClass(), e.getMessage());
                if (e.getCause() != null) {
                    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
                    final PrintStream ps = new PrintStream(baos);
                    e.getCause().printStackTrace(ps);
                    BaseTest.log(message);
                    BaseTest.log(baos.toString());
                }
                Assert.fail(message);
            }
        }

        if (expectedException != null && !gotAnException) {
            // we expected an exception but didn't see it
            Assert.fail(String.format("Test %s:%s expected exception %s but none was thrown", testClass.getClass().getSimpleName(), testName, expectedException.toString()));
        }
    }

    public static void assertMatchingFiles(final List<File> resultFiles, final List<String> expectedFiles, final boolean compareBamFilesSorted, final ValidationStringency stringency) throws IOException {
        Assert.assertEquals(resultFiles.size(), expectedFiles.size());
        for (int i = 0; i < resultFiles.size(); i++) {
            final File resultFile = resultFiles.get(i);
            final String expectedFileName = expectedFiles.get(i);
            final File expectedFile = new File(expectedFileName);
            if (expectedFileName.endsWith(".bam")){
                SamAssertionUtils.assertEqualBamFiles(resultFile, expectedFile, compareBamFilesSorted, stringency);
            } else {
                assertEqualTextFiles(resultFile, expectedFile);
            }
        }
    }

    public static void assertEqualTextFiles(final File resultFile, final File expectedFile) throws IOException {
        assertEqualTextFiles(resultFile, expectedFile, null);
    }

    /**
     * Compares two text files and ignores all lines that start with the comment prefix.
     */
    public static void assertEqualTextFiles(final File resultFile, final File expectedFile, final String commentPrefix) throws IOException {
        final Predicate<? super String> startsWithComment;
        if (commentPrefix == null){
            startsWithComment = s -> false;
        } else {
            startsWithComment = s -> s.startsWith(commentPrefix);
        }
        final List<String> actualLines = new XReadLines(resultFile).readLines().stream().filter(startsWithComment.negate()).collect(Collectors.toList());
        final List<String> expectedLines = new XReadLines(expectedFile).readLines().stream().filter(startsWithComment.negate()).collect(Collectors.toList());

        //For ease of debugging, we look at the lines first and only then check their counts
        int numUnequalLines = 0;
        final int minLen = Math.min(actualLines.size(), expectedLines.size());
        for (int i = 0; i < minLen; i++) {

            final String expectedLine = expectedLines.get(i);
            final String actualLine = actualLines.get(i);

            if ( !actualLine.equals(expectedLine) ) {
                logger.error( "Line number " + i + " (not counting comments) expected " +
                        expectedLine + " actual " + actualLine + '\n' +
                        "Expected :" + expectedLine  + '\n' +
                        "Actual   :" + actualLine  + '\n'
                );
                ++numUnequalLines;
            }
        }
        final boolean sizeMatches = actualLines.size() == expectedLines.size();

        // Check our error cases:
        if ( (numUnequalLines != 0) && (!sizeMatches) ) {
            throw new AssertionError("File sizes are unequal - actual = " + actualLines.size() + ", expected = " + expectedLines.size() + " AND detected unequal lines: " + numUnequalLines);
        }
        else if ( numUnequalLines != 0 ) {
            throw new AssertionError("Detected unequal lines: " + numUnequalLines);
        }
        else if (!sizeMatches) {
            throw new AssertionError("File sizes are unequal - actual = " + actualLines.size() + ", expected = " + expectedLines.size());
        }
    }

}
