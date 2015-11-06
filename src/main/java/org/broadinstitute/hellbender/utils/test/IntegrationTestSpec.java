package org.broadinstitute.hellbender.utils.test;

import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.picard.sam.SortSam;
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

    public IntegrationTestSpec(String args, List<String> expectedFileNames) {
        this.args = args;
        this.nOutputFiles = expectedFileNames.size();
        this.expectedException = null;
        this.expectedFileNames = expectedFileNames;
        this.compareBamFilesSorted = false;
    }

    public IntegrationTestSpec(String args, int nOutputFiles, Class<?> expectedException) {
        if (expectedException == null){
            throw new IllegalArgumentException("expected exception is null");
        }
        this.args = args;
        this.nOutputFiles = nOutputFiles;
        this.expectedException = expectedException;
        this.expectedFileNames = null;
        this.compareBamFilesSorted = false;
    }

    public boolean expectsException() {
        return expectedException != null;
    }

    public final Class<?> getExpectedException() {
        if (!expectsException())
            throw new GATKException("Tried to get exception for walker test that doesn't expect one");
        return expectedException;
    }

    public void setCompareBamFilesSorted(final boolean compareBamFilesSorted){
        this.compareBamFilesSorted = compareBamFilesSorted;
    }

    public String getArgs() {
        return args;
    }

    public Collection<String> expectedFileNames() {
        return expectedFileNames;
    }

    public void executeTest(final String name, CommandLineProgramTest test) throws IOException {
        List<File> tmpFiles = new ArrayList<>();
        for (int i = 0; i < nOutputFiles; i++) {
            String ext = DEFAULT_TEMP_EXTENSION;
            File fl = BaseTest.createTempFile(String.format(DEFAULT_TEMP_PREFIX + ".%d", i), ext);
            tmpFiles.add(fl);
        }

        final String args = String.format(getArgs(), tmpFiles.toArray());
        System.out.println(StringUtils.repeat('-', 80));

        if (expectsException()) {
            // this branch handles the case were we are testing that a walker will fail as expected
            executeTest(name, test, null, null, tmpFiles, args, getExpectedException());
        } else {
            List<String> expectedFileNames = new ArrayList<>();
            expectedFileNames.addAll(expectedFileNames());

            executeTest(name, test, null, expectedFileNames, tmpFiles, args, null);
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
    private void executeTest(String testName, CommandLineProgramTest testClass, File outputFileLocation, List<String> expectedFileNames, List<File> tmpFiles, String args, Class<?> expectedException) throws IOException {
        if (outputFileLocation != null) {
            args += " -O " + outputFileLocation.getAbsolutePath();
        }
        executeTest(testName, testClass, args, expectedException);

        if (expectedException == null && !expectedFileNames.isEmpty()) {
            assertMatchingFiles(tmpFiles, expectedFileNames, compareBamFilesSorted);
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
    private void executeTest(String testName, CommandLineProgramTest testClass, String args, Class<?> expectedException) {
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

    public static void assertMatchingFiles(final List<File> resultFiles, final List<String> expectedFiles, final boolean compareBamFilesSorted) throws IOException {
        Assert.assertEquals(resultFiles.size(), expectedFiles.size());
        for (int i = 0; i < resultFiles.size(); i++) {
            final File resultFile = resultFiles.get(i);
            final String expectedFileName = expectedFiles.get(i);
            final File expectedFile = new File(expectedFileName);
            if (expectedFileName.endsWith(".bam")){
                assertEqualBamFiles(resultFile, expectedFile, compareBamFilesSorted);
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

        Assert.assertEquals(actualLines.toString(), expectedLines.toString());
    }

    public static void assertEqualBamFiles(final File resultFile, final File expectedFile, final boolean compareBamFilesSorted) throws IOException {

        if (compareBamFilesSorted) {
            final File resultFileSorted= BaseTest.createTempFile("resultsFileSorted", ".bam");
            final File expectedFileSorted = BaseTest.createTempFile("expectedFileSorted", ".bam");

            sortSam(resultFile, resultFileSorted);
            sortSam(expectedFile, expectedFileSorted);

            SamAssertionUtils.assertSamsEqual(resultFileSorted, expectedFileSorted);
        } else {
            SamAssertionUtils.assertSamsEqual(resultFile, expectedFile);
        }
    }

    private static void sortSam(final File input, final File output) {
        final SortSam sort1 = new SortSam();
        sort1.INPUT = input;
        sort1.OUTPUT = output;
        sort1.SORT_ORDER = SAMFileHeader.SortOrder.coordinate;
        sort1.runTool();
    }

}
