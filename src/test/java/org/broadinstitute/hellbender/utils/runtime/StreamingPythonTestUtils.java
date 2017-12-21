package org.broadinstitute.hellbender.utils.runtime;

import org.testng.Assert;

public class StreamingPythonTestUtils {

    /**
     * Validate that we got a prompt written somewhere, either stdout out or stderr.
     */
    public static void assertPythonPrompt(final ProcessOutput po, final String prompt) {
        String stdOutOutput = null;
        String stdErrOutput = null;

        boolean bFoundPrompt = false;
        final StreamOutput stdOut = po.getStdout();
        if (stdOut != null) {
            stdOutOutput = stdOut.getBufferString();
            if (stdOutOutput.contains(prompt)) {
                bFoundPrompt = true;
            }
        }

        if (!bFoundPrompt) {
            final StreamOutput stderr = po.getStderr();
            if (stderr != null) {
                stdErrOutput = stderr.getBufferString();
                if (stdErrOutput.contains(prompt)) {
                    bFoundPrompt = true;
                }
            }
        }

        Assert.assertTrue(
                bFoundPrompt,
                String.format("Expected prompt [[%s]] but got stdout [[%s]] stderr [[%s]]",
                        prompt,
                        stdOutOutput == null ? "" : stdOutOutput,
                        stdErrOutput == null ? "" : stdErrOutput)
        );
    }

    /**
     * Validate that we got stdout output, not stderr. If exactMatch == false, then the output should start with
     * the expected value, otherwise it should match exactly.
     */
    public static void assertResponseOutput(final ProcessOutput po, final String expectedResult, final boolean requireExactMatch) {
        final StreamOutput stdOut = po.getStdout();
        if (stdOut == null) {
            final StreamOutput stderr = po.getStderr();
            if (stderr == null) {
                Assert.fail("Neither stdout nor stderr output was detected");
            } else {
                Assert.fail(
                        String.format(
                                "Expected stdout (%s), but got stderr: [[[%s]]]",
                                expectedResult,
                                stderr.getBufferString()));
            }
        } else {
            final String actualOutput = stdOut.getBufferString();
            boolean match = requireExactMatch ?
                    actualOutput.equals(expectedResult) :
                    actualOutput.contains(expectedResult);
            Assert.assertTrue(match, String.format("Expected \"[[[%s]]]\" but got \"[[[%s]]]\"", expectedResult, actualOutput));
        }
    }

}
