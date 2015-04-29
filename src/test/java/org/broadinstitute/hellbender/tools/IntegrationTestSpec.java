package org.broadinstitute.hellbender.tools;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.diffengine.DiffEngine;
import org.broadinstitute.hellbender.utils.read.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
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
import java.util.Collections;
import java.util.List;

public final class IntegrationTestSpec {
    public static final String DEFAULT_TEMP_EXTENSION = ".tmp";
    public static final String DEFAULT_TEMP_PREFIX = "walktest.tmp_param";

    private final String args;
    private final Class<?> expectedException;
    private final int nOutputFiles;
    private final List<String> expectedFileNames;

    public IntegrationTestSpec(final String args, final List<String> expectedFileNames) {
        this.args = args;
        this.nOutputFiles = expectedFileNames.size();
        this.expectedException = null;
        this.expectedFileNames = expectedFileNames;
    }

    public IntegrationTestSpec(final String args, final int nOutputFiles, final Class<?> expectedException) {
        if (expectedException == null){
            throw new IllegalArgumentException("expected exception is null");
        }
        this.args = args;
        this.nOutputFiles = nOutputFiles;
        this.expectedException = expectedException;
        this.expectedFileNames = null;
    }

    public boolean expectsException() {
        return expectedException != null;
    }

    public final Class<?> getExpectedException() {
        if (!expectsException()) {
            throw new GATKException("Tried to get exception for walker test that doesn't expect one");
        }
        return expectedException;
    }

    public String getArgs() {
        return args;
    }

    public Collection<String> expectedFileNames() {
        return expectedFileNames;
    }

    public void executeTest(final String name, final CommandLineProgramTest test) throws IOException {
        final List<File> tmpFiles = new ArrayList<>();
        for (int i = 0; i < nOutputFiles; i++) {
            final String ext = DEFAULT_TEMP_EXTENSION;
            final File fl = BaseTest.createTempFile(String.format(DEFAULT_TEMP_PREFIX + ".%d", i), ext);
            tmpFiles.add(fl);
        }

        final String args = String.format(getArgs(), tmpFiles.toArray());
        System.out.println(StringUtils.repeat('-', 80));

        if (expectsException()) {
            // this branch handles the case were we are testing that a walker will fail as expected
            executeTest(name, test, null, null, tmpFiles, args, getExpectedException());
        } else {
            final List<String> expectedFileNames = new ArrayList<>();
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
     * @return a pair of file and string lists
     */
    private void executeTest(final String testName, final CommandLineProgramTest testClass, final File outputFileLocation, final List<String> expectedFileNames, final List<File> tmpFiles, String args, final Class<?> expectedException) throws IOException {
        if (outputFileLocation != null) {
            args += " -O " + outputFileLocation.getAbsolutePath();
        }
        executeTest(testName, testClass, args, expectedException);

        if (expectedException == null && !expectedFileNames.isEmpty()) {
            assertMatchingFiles(tmpFiles, expectedFileNames);
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
    private void executeTest(final String testName, final CommandLineProgramTest testClass, final String args, final Class<?> expectedException) {
        final String[] command = Utils.escapeExpressions(args);
        // run the executable
        boolean gotAnException = false;
        try {
            final String now = LocalDateTime.now().format(DateTimeFormatter.ofPattern("HH:mm:ss"));
            final String cmdline = String.join(" ", command);
            System.out.println(String.format("[%s] Executing test %s:%s with GATK arguments: %s", now, testClass.getClass().getSimpleName(), testName, cmdline));
            // also write the command line to the HTML log for convenient follow-up
            // do the replaceAll so paths become relative to the current
            BaseTest.log(cmdline.replaceAll(BaseTest.publicTestDirRoot, ""));
            final Object result = testClass.runCommandLine(command);
        } catch (final Exception e) {
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

    public static void assertMatchingFiles(final List<File> resultFiles, final List<String> expectedFiles) throws IOException {
        Assert.assertEquals(resultFiles.size(), expectedFiles.size());
        for (int i = 0; i < resultFiles.size(); i++) {
            final File resultFile = resultFiles.get(i);
            final String expectedFileName = expectedFiles.get(i);
            final File expectedFile = new File(expectedFileName);
            if (expectedFileName.endsWith(".bam")){
                compareBamFiles(resultFile, expectedFile);
            } else if (expectedFileName.endsWith(".vcf")){
                compareVcfFiles(resultFile, expectedFile);
            } else {
                compareTextFiles(resultFile, expectedFile);
            }
        }
    }

    private static void compareVcfFiles(final File resultFile, final File expectedFile) throws IOException {
        final List<String> actualLines = new XReadLines(resultFile).readLines();
        final List<String> expectedLines = new XReadLines(expectedFile).readLines();
        if (!actualLines.equals(expectedLines)){
             new DiffEngine().diffFiles(resultFile, expectedFile);
        }
    }

    public static void compareTextFiles(final File resultFile, final File expectedFile) throws IOException {
        final List<String> actualLines = new XReadLines(resultFile).readLines();
        final List<String> expectedLines = new XReadLines(expectedFile).readLines();
        Assert.assertEquals(actualLines.toString(), expectedLines.toString());
    }

    public static void compareBamFiles(final File resultFile, final File expectedFile) throws IOException {
        final boolean equal = SamAssertionUtils.areSamsEqual(resultFile, expectedFile);
        if (!equal){
            new DiffEngine().diffFiles(resultFile, expectedFile);
        }
    }
}
