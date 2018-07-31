package org.broadinstitute.hellbender.utils.runtime;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

public final class ProcessControllerUnitTest extends GATKBaseTest {
    private static final String NL = String.format("%n");

    @Test
    public void testReuseAfterError() {
        ProcessController controller = new ProcessController();

        ProcessSettings job;

        for (int i = 0; i < 3; i++) {
            // Test bad command
            job = new ProcessSettings(new String[] {"no_such_command"});
            try {
                controller.exec(job);
            } catch (GATKException e) {
                /* Was supposed to throw an exception */
            }

            // Test exit != 0
            job = new ProcessSettings(new String[] {"cat", "non_existent_file"});
            int exitValue = controller.exec(job).getExitValue();
            Assert.assertTrue(exitValue != 0, "'cat' non existent file returned 0");

            // Text success
            job = new ProcessSettings(new String[] {"echo", "Hello World"});
            exitValue = controller.exec(job).getExitValue();
            Assert.assertEquals(exitValue, 0, "Echo failed");
        }
    }

    @Test
    public void testEnvironment() {
        String key = "MY_NEW_VAR";
        String value = "value is here";

        ProcessSettings job = new ProcessSettings(new String[] {"sh", "-c", "echo $"+key});
        job.getStdoutSettings().setBufferSize(-1);
        job.setRedirectErrorStream(true);

        Map<String, String> env = new LinkedHashMap<>(System.getenv());
        env.put(key, value);
        job.setEnvironment(env);

        ProcessController controller = new ProcessController();
        ProcessOutput result = controller.exec(job);
        int exitValue = result.getExitValue();

        Assert.assertEquals(exitValue, 0, "Echo environment variable failed");
        Assert.assertEquals(result.getStdout().getBufferString(), value + NL, "Echo environment returned unexpected output");
    }

    @Test
    public void testDirectory() throws IOException {
        File dir = createTempDir("temp.").getCanonicalFile();

        ProcessSettings job = new ProcessSettings(new String[] {"pwd"});
        job.getStdoutSettings().setBufferSize(-1);
        job.setRedirectErrorStream(true);
        job.setDirectory(dir);

        ProcessController controller = new ProcessController();
        ProcessOutput result = controller.exec(job);
        int exitValue = result.getExitValue();

        Assert.assertEquals(exitValue, 0, "Getting working directory failed");

        Assert.assertEquals(result.getStdout().getBufferString(), dir.getAbsolutePath() + NL,
                "Setting/getting working directory returned unexpected output");
    }

    @Test
    public void testReadStdInBuffer() {
        String bufferText = "Hello from buffer";
        ProcessSettings job = new ProcessSettings(new String[] {"cat"});
        job.getStdoutSettings().setBufferSize(-1);
        job.setRedirectErrorStream(true);
        job.getStdinSettings().setInputBuffer(bufferText);

        ProcessController controller = new ProcessController();
        ProcessOutput output = controller.exec(job);

        Assert.assertEquals(output.getStdout().getBufferString(), bufferText,
                "Unexpected output from cat stdin buffer");
    }

    @Test
    public void testReadStdInFile() {
        File input = null;
        try {
            String fileText = "Hello from file";
            input = IOUtils.writeTempFile(fileText, "stdin.", ".txt");

            ProcessSettings job = new ProcessSettings(new String[] {"cat"});
            job.getStdoutSettings().setBufferSize(-1);
            job.setRedirectErrorStream(true);
            job.getStdinSettings().setInputFile(input);

            ProcessController controller = new ProcessController();
            ProcessOutput output = controller.exec(job);

            Assert.assertEquals(output.getStdout().getBufferString(), fileText,
                    "Unexpected output from cat stdin file");
        } finally {
            FileUtils.deleteQuietly(input);
        }
    }

    @Test
    public void testWriteStdOut() {
        final String testInput = "Testing to stdout";
        final ProcessSettings job = new ProcessSettings(new String[] {"echo", testInput});
        job.getStdoutSettings().printStandard(true);
        job.getStdoutSettings().setBufferSize(-1);
        job.setRedirectErrorStream(true);

        final ProcessController controller = new ProcessController();
        ProcessOutput processOutput = controller.exec(job);
        Assert.assertEquals(processOutput.getExitValue(), 0);
        Assert.assertEquals(processOutput.getStdout().getBufferString(), testInput + NL);
    }

    @Test
    public void testErrorToOut() throws IOException {
        File outFile = null;
        File errFile = null;
        try {
            outFile = GATKBaseTest.createTempFile("temp", "");
            errFile = GATKBaseTest.createTempFile("temp", "");

            ProcessSettings job = new ProcessSettings(new String[]{"cat", "non_existent_file"});
            job.getStdoutSettings().setOutputFile(outFile);
            job.getStdoutSettings().setBufferSize(-1);
            job.getStderrSettings().setOutputFile(errFile);
            job.getStderrSettings().setBufferSize(-1);
            job.setRedirectErrorStream(true);

            ProcessOutput result = new ProcessController().exec(job);
            int exitValue = result.getExitValue();

            Assert.assertTrue(exitValue != 0, "'cat' non existent file returned 0");

            String fileString, bufferString;

            fileString = FileUtils.readFileToString(outFile, StandardCharsets.UTF_8);
            Assert.assertTrue(fileString.length() > 0, "Out file was length 0");

            bufferString = result.getStdout().getBufferString();
            Assert.assertTrue(bufferString.length() > 0, "Out buffer was length 0");

            Assert.assertFalse(result.getStdout().isBufferTruncated(), "Out buffer was truncated");
            Assert.assertEquals(bufferString.length(), fileString.length(), "Out buffer length did not match file length");

            fileString = FileUtils.readFileToString(errFile, StandardCharsets.UTF_8);
            Assert.assertEquals(fileString, "", "Unexpected output to err file");

            bufferString = result.getStderr().getBufferString();
            Assert.assertEquals(bufferString, "", "Unexepected output to err buffer");
        } finally {
            FileUtils.deleteQuietly(outFile);
            FileUtils.deleteQuietly(errFile);
        }
    }

    @Test
    public void testErrorToErr() throws IOException {
        File outFile = null;
        File errFile = null;
        try {
            outFile = GATKBaseTest.createTempFile("temp", "");
            errFile = GATKBaseTest.createTempFile("temp", "");

            ProcessSettings job = new ProcessSettings(new String[]{"cat", "non_existent_file"});
            job.getStdoutSettings().setOutputFile(outFile);
            job.getStdoutSettings().setBufferSize(-1);
            job.getStderrSettings().setOutputFile(errFile);
            job.getStderrSettings().setBufferSize(-1);
            job.setRedirectErrorStream(false);

            ProcessOutput result = new ProcessController().exec(job);
            int exitValue = result.getExitValue();

            Assert.assertTrue(exitValue != 0, "'cat' non existent file returned 0");

            String fileString, bufferString;

            fileString = FileUtils.readFileToString(errFile, StandardCharsets.UTF_8);
            Assert.assertTrue(fileString.length() > 0, "Err file was length 0");

            bufferString = result.getStderr().getBufferString();
            Assert.assertTrue(bufferString.length() > 0, "Err buffer was length 0");

            Assert.assertFalse(result.getStderr().isBufferTruncated(), "Err buffer was truncated");
            Assert.assertEquals(bufferString.length(), fileString.length(), "Err buffer length did not match file length");

            fileString = FileUtils.readFileToString(outFile, StandardCharsets.UTF_8);
            Assert.assertEquals(fileString, "", "Unexpected output to out file");

            bufferString = result.getStdout().getBufferString();
            Assert.assertEquals(bufferString, "", "Unexepected output to out buffer");
        } finally {
            FileUtils.deleteQuietly(outFile);
            FileUtils.deleteQuietly(errFile);
        }
    }

    private static final String TRUNCATE_TEXT = "Hello World";
    private static final byte[] TRUNCATE_OUTPUT_BYTES = (TRUNCATE_TEXT + NL).getBytes();

    /**
     * @return Test truncating content vs. not truncating (run at -1/+1 size)
     */
    @DataProvider(name = "truncateSizes")
    public Object[][] getTruncateBufferSizes() {
        int l = TRUNCATE_OUTPUT_BYTES.length;
        return new Object[][]{
                new Object[]{0, 0},
                new Object[]{l, l},
                new Object[]{l + 1, l},
                new Object[]{l - 1, l - 1}
        };
    }

    @Test(dataProvider = "truncateSizes")
    public void testTruncateBuffer(int truncateLen, int expectedLen) {
        byte[] expected = Arrays.copyOf(TRUNCATE_OUTPUT_BYTES, expectedLen);

        String[] command = {"echo", TRUNCATE_TEXT};
        ProcessController controller = new ProcessController();

        ProcessSettings job = new ProcessSettings(command);
        job.getStdoutSettings().setBufferSize(truncateLen);
        ProcessOutput result = controller.exec(job);

        int exitValue = result.getExitValue();

        Assert.assertEquals(exitValue, 0,
                String.format("Echo returned %d: %s", exitValue, TRUNCATE_TEXT));

        byte[] bufferBytes = result.getStdout().getBufferBytes();

        Assert.assertEquals(bufferBytes, expected,
                String.format("Output buffer didn't match (%d vs %d)", expected.length, bufferBytes.length));

        boolean truncated = result.getStdout().isBufferTruncated();

        Assert.assertEquals(truncated, TRUNCATE_OUTPUT_BYTES.length > truncateLen,
                "Unexpected buffer truncation result");
    }

    private static final String[] LONG_COMMAND = getLongCommand();
    private static final String LONG_COMMAND_STRING = StringUtils.join(LONG_COMMAND, " ");
    private static final String LONG_COMMAND_DESCRIPTION = "<long command>";

    @DataProvider(name = "echoCommands")
    public Object[][] getEchoCommands() {

        new EchoCommand(new String[]{"echo", "Hello", "World"}, "Hello World" + NL);
        new EchoCommand(new String[]{"echo", "'Hello", "World"}, "'Hello World" + NL);
        new EchoCommand(new String[]{"echo", "Hello", "World'"}, "Hello World'" + NL);
        new EchoCommand(new String[]{"echo", "'Hello", "World'"}, "'Hello World'" + NL);

        String[] longCommand = new String[LONG_COMMAND.length + 1];
        longCommand[0] = "echo";
        System.arraycopy(LONG_COMMAND, 0, longCommand, 1, LONG_COMMAND.length);
        new EchoCommand(longCommand, LONG_COMMAND_STRING + NL) {
            @Override
            public String toString() {
                return LONG_COMMAND_DESCRIPTION;
            }
        };

        return TestDataProvider.getTests(EchoCommand.class);
    }

    @Test(dataProvider = "echoCommands")
    public void testEcho(EchoCommand script) throws IOException {
        File outputFile = null;
        try {
            outputFile = GATKBaseTest.createTempFile("temp", "");

            ProcessSettings job = new ProcessSettings(script.command);
            if (script.output != null) {
                job.getStdoutSettings().setOutputFile(outputFile);
                job.getStdoutSettings().setBufferSize(script.output.getBytes().length);
            }

            ProcessOutput result = new ProcessController().exec(job);
            int exitValue = result.getExitValue();

            Assert.assertEquals(exitValue, 0,
                    String.format("Echo returned %d: %s", exitValue, script));

            if (script.output != null) {

                String fileString = FileUtils.readFileToString(outputFile, StandardCharsets.UTF_8);
                Assert.assertEquals(fileString, script.output,
                        String.format("Output file didn't match (%d vs %d): %s",
                                fileString.length(), script.output.length(), script));

                String bufferString = result.getStdout().getBufferString();
                Assert.assertEquals(bufferString, script.output,
                        String.format("Output content didn't match (%d vs %d): %s",
                                bufferString.length(), script.output.length(), script));

                Assert.assertFalse(result.getStdout().isBufferTruncated(),
                        "Output content was truncated: " + script);
            }
        } finally {
            FileUtils.deleteQuietly(outputFile);
        }
    }

    @Test(expectedExceptions = GATKException.class)
    public void testUnableToStart() {
        ProcessSettings job = new ProcessSettings(new String[]{"no_such_command"});
        new ProcessController().exec(job);
    }

    @DataProvider(name = "scriptCommands")
    public Object[][] getScriptCommands() {
        new ScriptCommand(true, "echo Hello World", "Hello World" + NL);
        new ScriptCommand(false, "echo 'Hello World", null);
        new ScriptCommand(false, "echo Hello World'", null);
        new ScriptCommand(true, "echo 'Hello World'", "Hello World" + NL);
        new ScriptCommand(true, "echo \"Hello World\"", "Hello World" + NL);
        new ScriptCommand(false, "no_such_echo Hello World", null);
        new ScriptCommand(true, "echo #", NL);
        new ScriptCommand(true, "echo \\#", "#" + NL);
        new ScriptCommand(true, "echo \\\\#", "\\#" + NL);

        new ScriptCommand(true, "echo " + LONG_COMMAND_STRING, LONG_COMMAND_STRING + NL) {
            @Override
            public String toString() {
                return LONG_COMMAND_DESCRIPTION;
            }
        };

        return TestDataProvider.getTests(ScriptCommand.class);
    }

    @Test(dataProvider = "scriptCommands")
    public void testScript(ScriptCommand script) throws IOException {
        File scriptFile = null;
        File outputFile = null;
        try {
            scriptFile = writeScript(script.content);
            outputFile = GATKBaseTest.createTempFile("temp", "");

            ProcessSettings job = new ProcessSettings(new String[]{"sh", scriptFile.getAbsolutePath()});
            if (script.output != null) {
                job.getStdoutSettings().setOutputFile(outputFile);
                job.getStdoutSettings().setBufferSize(script.output.getBytes().length);
            }

            ProcessOutput result = new ProcessController().exec(job);
            int exitValue = result.getExitValue();

            Assert.assertEquals(exitValue == 0, script.succeed,
                    String.format("Script returned %d: %s", exitValue, script));

            if (script.output != null) {

                String fileString = FileUtils.readFileToString(outputFile, StandardCharsets.UTF_8);
                Assert.assertEquals(fileString, script.output,
                        String.format("Output file didn't match (%d vs %d): %s",
                                fileString.length(), script.output.length(), script));

                String bufferString = result.getStdout().getBufferString();
                Assert.assertEquals(bufferString, script.output,
                        String.format("Output content didn't match (%d vs %d): %s",
                                bufferString.length(), script.output.length(), script));

                Assert.assertFalse(result.getStdout().isBufferTruncated(),
                        "Output content was truncated: " + script);
            }
        } finally {
            FileUtils.deleteQuietly(scriptFile);
            FileUtils.deleteQuietly(outputFile);
        }
    }

    private static String[] getLongCommand() {
        // This command fails on some systems with a 4096 character limit when run via the old sh -c "echo ...",
        // but works on the same systems when run via sh <script>
        int cnt = 500;
        String[] command = new String[cnt];
        for (int i = 1; i <= cnt; i++) {
            command[i - 1] = String.format("%03d______", i);
        }
        return command;
    }

    private static File writeScript(String contents) {
        try {
            File file = GATKBaseTest.createTempFile("temp", "");
            FileUtils.writeStringToFile(file, contents, StandardCharsets.UTF_8);
            return file;
        } catch (IOException e) {
            throw new UserException.BadTempDir(e.getMessage(), e);
        }
    }

    private static class EchoCommand extends TestDataProvider {
        public final String[] command;
        public final String output;

        public EchoCommand(String[] command, String output) {
            super(EchoCommand.class);
            this.command = command;
            this.output = output;
        }

        @Override
        public String toString() {
            return StringUtils.join(command, " ");
        }
    }

    public static class ScriptCommand extends TestDataProvider {
        public final boolean succeed;
        public final String content;
        public final String output;

        public ScriptCommand(boolean succeed, String content, String output) {
            super(ScriptCommand.class);
            this.succeed = succeed;
            this.content = content;
            this.output = output;
        }

        @Override
        public String toString() {
            return content;
        }
    }
}
