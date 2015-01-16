package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.*;

public final class WalkerTestSpec {
    private final String args;
    private final Class<?> expectedException;
    private final int nOutputFiles;
    private final List<String> expectedFileNames;

    List<String> exts = null;

    // the default output path for the integration test
    private File outputFileLocation = null;


    public WalkerTestSpec(String args, List<String> expectedFileNames) {
        this.args = args;
        this.nOutputFiles = expectedFileNames.size();
        this.expectedException = null;
        this.expectedFileNames = expectedFileNames;
    }

    public WalkerTestSpec(String args, int nOutputFiles, Class<?> expectedException) {
        this.args = args;
        this.nOutputFiles = nOutputFiles;
        this.expectedException = expectedException;
        this.expectedFileNames = null;
    }

    public boolean expectsException() {
        return expectedException != null;
    }

    public Class<?> getExpectedException() {
        if (!expectsException())
            throw new GATKException("Tried to get exception for walker test that doesn't expect one");
        return expectedException;
    }

    public String getArgs() {
        return args;
    }

    public void setOutputFileLocation(File outputFileLocation) {
        this.outputFileLocation = outputFileLocation;
    }

    protected File getOutputFileLocation() {
        return outputFileLocation;
    }

    public Collection<String> expectedFileNames() {
        return expectedFileNames;
    }

    public void executeTest(final String name, CommandLineProgramTest test) throws IOException {
        List<File> tmpFiles = new ArrayList<>();
        for (int i = 0; i < nOutputFiles; i++) {
            String ext = exts == null ? ".tmp" : "." + exts.get(i);
            File fl = BaseTest.createTempFile(String.format("walktest.tmp_param.%d", i), ext);
            tmpFiles.add(fl);
        }

        final String args = String.format(getArgs(), tmpFiles.toArray());
        System.out.println(Utils.dupString('-', 80));

        if (expectsException()) {
            // this branch handles the case were we are testing that a walker will fail as expected
            executeTest(name, test, getOutputFileLocation(), null, tmpFiles, args, getExpectedException());
        } else {
            List<String> expectedFileNames = new LinkedList<>();
            expectedFileNames.addAll(expectedFileNames());

            executeTest(name, test, getOutputFileLocation(), expectedFileNames, tmpFiles, args, null);
        }
    }

    /**
     * execute the test, given the following:
     *
     * @param testName              the name of the test
     * @param testClass             the object that contains the test
     * @param expectedTextFileNames the list of expectedFileNames
     * @param tmpFiles              the temp file corresponding to the expectedFileNames list
     * @param args                  the argument list
     * @param expectedException     the expected exception or null
     * @return a pair of file and string lists
     */
    private void executeTest(String testName, CommandLineProgramTest testClass, File outputFileLocation, List<String> expectedTextFileNames, List<File> tmpFiles, String args, Class<?> expectedException) throws IOException {
        if (outputFileLocation != null)
            args += " -o " + outputFileLocation.getAbsolutePath();
        executeTest(testName, testClass, args, expectedException);

        if (expectedException == null) {
            assertMatchingFiles(tmpFiles, expectedTextFileNames);
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
            final String now = new SimpleDateFormat("HH:mm:ss").format(new Date());
            final String cmdline = Utils.join(" ", command);
            System.out.println(String.format("[%s] Executing test %s:%s with GATK arguments: %s", now, testClass.getClass().getSimpleName(), testName, cmdline));
            // also write the command line to the HTML log for convenient follow-up
            // do the replaceAll so paths become relative to the current
            BaseTest.log(cmdline.replaceAll(BaseTest.publicTestDirRoot, ""));
            Object result = testClass.runCommandLine(command);
        } catch (Exception e) {
            gotAnException = true;
            if (expectedException != null) {
                // we expect an exception
                //System.out.println(String.format("Wanted exception %s, saw %s", expectedException, e.getClass()));
                if (expectedException.isInstance(e)) {
                    // it's the type we expected
                    //System.out.println(String.format("  => %s PASSED", name));
                } else {
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
            } else {
                // we didn't expect an exception but we got one :-(
                throw new RuntimeException(e);
            }
        }

        if (expectedException != null && !gotAnException) {
            // we expected an exception but didn't see it
            Assert.fail(String.format("Test %s:%s expected exception %s but none was thrown", testClass.getClass().getSimpleName(), testName, expectedException.toString()));
        }
    }

    public void assertMatchingFiles(List<File> resultFiles, List<String> expectedFiles) throws IOException {
        Assert.assertEquals(resultFiles.size(), expectedFiles.size());
        for (int i = 0; i < resultFiles.size(); i++) {
            File resultFile = resultFiles.get(i);
            String expectedFileName = expectedFiles.get(i);
            File expectedFile = new File(expectedFileName);
            if (expectedFileName.endsWith(".bam")){
                compareBamFiles(resultFile, expectedFile);
            } else {
                compareTextFiles(resultFile, expectedFile);
            }
        }
    }

    private void compareTextFiles(File resultFile, File expectedFile) throws IOException {
        List<String> actualLines = new XReadLines(resultFile).readLines();
        List<String> expectedLines = new XReadLines(expectedFile).readLines();
        Assert.assertEquals(actualLines.toString(), expectedLines.toString());
    }

    private void compareBamFiles(File resultFile, File expectedFile) throws IOException {
//     //TODO - there's a bug in CompareSAMs that prevents us from using it.
       // Until it's fixed, we're converting bams to sams and comparing strings.
//     final CompareSAMs compareSAMs = new CompareSAMs();
//     compareSAMs.samFiles = Arrays.asList(resultFile, expectedFile);
//     compareSAMs.doWork();
//     Assert.assertTrue(compareSAMs.areEqual(), "actual:" + resultFile.getAbsolutePath() + " expected:" + expectedFile);

        //HACK comparing bam files as text files, first converting them to sam files
        PrintReads pr = new PrintReads();
        File expectedSAMfile = BaseTest.createTempFile("expected", ".sam");
        pr.OUTPUT=  expectedSAMfile;
        pr.READS_FILES = Arrays.asList(expectedFile);
        pr.runTool();

        PrintReads pr1 = new PrintReads();
        File actualSAMfile = BaseTest.createTempFile("actual", ".sam");
        pr1.OUTPUT= actualSAMfile;
        pr1.READS_FILES = Arrays.asList(resultFile);
        pr1.runTool();

        compareTextFiles(actualSAMfile, expectedSAMfile);
    }
}
