package org.broadinstitute.hellbender.utils.python;

import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.testng.TestException;

import java.util.List;

/**
 * A test runner for Python modules that implement the unittest testing framework
 */
public class PythonUnitTestRunner extends StreamingPythonScriptExecutor<String> {
    private static final String NL = System.lineSeparator();
    private static final String TAB = "\t";

    /**
     * The start method must be called to actually start the remote executable.
     */
    public PythonUnitTestRunner() {
        super(true);
    }

    /**
     * Initialize the test runner
     * @param imports any python packages necessary for the test, usually the package the test belongs to
     * @param additionalCommands other pre-test python commands
     * @param pythonProcessArgs arguments to use when invoking the python interpreter
     */
    public void setup(final List<String> imports, final List<String> additionalCommands, final List<String> pythonProcessArgs) {
        start(pythonProcessArgs);
        for (final String packageName : imports) {
            sendSynchronousCommand("import " + packageName + System.lineSeparator());
        }
        for (final String command : additionalCommands) {
            sendSynchronousCommand(command + System.lineSeparator());
        }
    }

    /**
     * Java method to run a python unittest and report results to stdout, throwing a TestNG TestException if applicable
     * @param testModuleName module name relative to src/main/python/org/broadinstitute/hellbender/
     */
    public void runTest(final String testModuleName) {
        try {
            sendSynchronousCommand("from unittest import main" + System.lineSeparator());
            sendSynchronousCommand("test = main(module=" + testModuleName + ", exit=False)"  + System.lineSeparator());
            printResults(TestResultType.SKIPPED);
            printResults(TestResultType.ERRORS);
            printResults(TestResultType.FAILURES);
            printResults(TestResultType.UNEXPECTED_SUCCESSES);  //TODO: this seems to output null if there were no tests with expected exceptions
            passOrFail();
        } catch (TestException e) {
            throw e;
        } finally {
            terminate();
        }
    }

    /**
     * Actually throw an exception if tests errored or failed; to be called after printResults
     */
    private void passOrFail() {
        try {
            sendSynchronousCommand("assert(len(test.result.errors) == 0)" + System.lineSeparator());
        } catch (final PythonScriptExecutorException e) {
            throw new TestException("Test errors occurred.  See stderr for details.");
        }

        try {
            sendSynchronousCommand("assert(len(test.result.failures) == 0)" + System.lineSeparator());
        } catch (final PythonScriptExecutorException e) {
            throw new TestException("Test failures occurred.  See stderr for details.");
        }

        try {
            sendSynchronousCommand("assert(len(test.result.unexpectedSuccesses) == 0)" + System.lineSeparator());
        } catch (final PythonScriptExecutorException e) {
            throw new TestException("Unexpected successes occurred.  See stderr for details.");
        }
    }

    /**
     *  Report back the count of issues of this type, the relevant test names, and messages
     * @param type type of issue encountered when running tests
     */
    private void printResults(final TestResultType type) {
        ProcessOutput out;
        out = sendSynchronousCommand("print(str(len(test.result." + type.getFieldName() + ")) + \" " + type.getDescription() +"\")" + NL);
        System.err.println(out.getStdout());
        final StringBuilder sb = new StringBuilder();
        sb.append("if len(test.result." + type.getFieldName() + ") > 0:" + NL);
        sb.append(TAB + "for result_pair in test.result." + type.getFieldName() + ":" + NL);
        sb.append(TAB + TAB + "print(result_pair[0].id())" + NL );
        sb.append(TAB + TAB + "print(result_pair[1])" + NL + NL + NL); // extra newlines for if and for loop
        try {
            out = sendSynchronousCommand(sb.toString());
            if (out.getStdout() != null) {
                System.err.println(out.getStdout());
            }
        } catch (final PythonScriptExecutorException e) {
            System.err.println("Test result retrieval failed: " + e);
        }
    }

    /**
     * Classifications of test result types within the python unittest framework
     * Additional types exist: https://docs.python.org/3/library/unittest.html#unittest.TestResult
     */
    private enum TestResultType {
        ERRORS("errors", "test errors"),
        FAILURES("failures", "test failures"),
        SKIPPED("skipped", "skipped tests"),
        EXPECTED_FAILURES("expectedFailures", "expected failures"),
        UNEXPECTED_SUCCESSES("unexpectedSuccesses", "unexpected successes");

        private final String fieldName;
        private final String description;

        TestResultType(final String pythonField, final String description) {
            fieldName = pythonField;
            this.description = description;
        }

       public String getFieldName() {
            return fieldName;
        }

       public String getDescription() {return description;}
    }
}