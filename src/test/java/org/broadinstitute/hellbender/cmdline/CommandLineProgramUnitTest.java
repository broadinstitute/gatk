package org.broadinstitute.hellbender.cmdline;


import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class CommandLineProgramUnitTest extends GATKBaseTest {

    @Test
    public void testGetUsage(){
        final CommandLineProgram clp = getClp();
        String usage = clp.getUsage();
        GATKBaseTest.assertContains(usage, "Usage:");
    }

    @Test
    public void testGetCommandLine(){
        final CommandLineProgram clp = getClp();
        Assert.assertNull(clp.getCommandLine()); //should be null since no args were specified
        Assert.assertFalse(clp.parseArgs(new String[]{"--" + SpecialArgumentsCollection.HELP_FULLNAME}));
        assertContains(clp.getCommandLine(), SpecialArgumentsCollection.HELP_FULLNAME); //now it should be filled in
    }

    private static class ValidationFailer extends CommandLineProgram{
        public static final String ERROR1 = "first error";
        public static final String ERROR2 = "second error";

        @Override
        protected Object doWork() {
            return null;
        }

        @Override
        protected String[] customCommandLineValidation(){
            return new String[]{ERROR1, ERROR2};
        }
    }

    @Test
    public void testCustomValidationFailThrowsCommandLineException(){
        ValidationFailer clp = new ValidationFailer();
        try{
            clp.parseArgs(new String[] {});
            Assert.fail("Should have thrown an exception");
        } catch (final CommandLineException e){
            final String message = e.getMessage();
            assertContains(message, ValidationFailer.ERROR1);
            assertContains(message, ValidationFailer.ERROR2);
        }
    }

    private static CommandLineProgram getClp() {
        return new CommandLineProgram() {
            @Override
            protected Object doWork() {
                return null;
            }
        };
    }

    @CommandLineProgramProperties(
            summary = "TestGATKToolWithSuppressedFileExpansion",
            oneLineSummary = "TestGATKToolWithSuppressedFileExpansion",
            programGroup = TestProgramGroup.class
    )
    private static final class TestGATKToolWithSuppressedFileExpansion extends CommandLineProgram {

        @Argument(fullName = "suppressedExpansionArg", suppressFileExpansion = true, optional = true)
        List<String> suppressedExpansionArg = new ArrayList<>();

        @Argument(fullName = "notSuppressedExpansionArg", optional = true)
        List<String> notSuppressedExpansionArg = new ArrayList<>();

        @Override
        public Object doWork() {
            //no op
            return 1;
        }

        public List<String> getSuppressedExpansionArg() {
            return suppressedExpansionArg;
        }

        public List<String> getNotSuppressedExpansionArg() {
            return notSuppressedExpansionArg;
        }
    }

    @DataProvider(name = "commandLineExpansionExtensions")
    private Object[][] commandLineExpansionExtensions() throws IOException {
        final String listFileContent = "this is a .list file";
        final String argsFileContent = "this is a .args file";
        return new Object[][]{
                {createTempExpansionTestFile(".list", listFileContent), listFileContent},
                {createTempExpansionTestFile(".args", argsFileContent), argsFileContent},
        };
    }

    @Test(dataProvider = "commandLineExpansionExtensions")
    public void testListExpansionSuppression(
            final File testFile,
            @SuppressWarnings("unused") String fileContents) throws IOException {
        final TestGATKToolWithSuppressedFileExpansion tool = new TestGATKToolWithSuppressedFileExpansion();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final String[] args = {
                "--suppressedExpansionArg", testFile.getCanonicalPath()
        };
        clp.parseArguments(System.out, args);

        // for the arg with file expansion suppressed, the list should contain the actual file name, not the contents
        Assert.assertEquals(
                tool.getSuppressedExpansionArg(),
                Collections.singletonList(testFile.getCanonicalPath()));
    }

    @Test(dataProvider = "commandLineExpansionExtensions")
    public void testListExpansion(final File testFile, final String fileContents) throws IOException {
        final TestGATKToolWithSuppressedFileExpansion tool = new TestGATKToolWithSuppressedFileExpansion();
        final CommandLineParser clp = new CommandLineArgumentParser(tool);
        final String[] args = {
                "--notSuppressedExpansionArg", testFile.getCanonicalPath()
        };
        clp.parseArguments(System.out, args);

        // for the arg with file expansion NOT suppressed (defualt case), the list should contain the file contents
        Assert.assertEquals(
                tool.getNotSuppressedExpansionArg(),
                Collections.singletonList(fileContents));
    }

    private File createTempExpansionTestFile(final String extension, final String fileContentLine) throws IOException {
        File tempFile = createTempFile("testExpansionSuppression", extension);
        try (final FileWriter fileWriter = new FileWriter(tempFile)) {
            fileWriter.write(fileContentLine);
        }
        return tempFile;
    }

    @CommandLineProgramProperties(
            summary = "ExperimentalTool",
            oneLineSummary = "ExperimentalTool",
            programGroup = TestProgramGroup.class
    )
    @ExperimentalFeature
    private static class TestExperimentalTool extends CommandLineProgram {

        @Override
        protected Object doWork() {
            return null;
        }
    }

    @CommandLineProgramProperties(
            summary = "BetaTool",
            oneLineSummary = "BetaTool",
            programGroup = TestProgramGroup.class
    )
    @BetaFeature
    private static class TestBetaTool extends CommandLineProgram {
        @Override
        protected Object doWork() {
            return null;
        }
    }

    @CommandLineProgramProperties(
            summary = "ProductionTool",
            oneLineSummary = "ProductionTool",
            programGroup = TestProgramGroup.class
    )
    private static class TestProductionTool extends CommandLineProgram {
        @Override
        protected Object doWork() {
            return null;
        }
    }

    @DataProvider(name = "toolMaturityLevels")
    public Object[][] getToolMaturityLevels() {
        return new Object[][]{
                {new TestExperimentalTool(), true, false, "EXPERIMENTAL"},
                {new TestBetaTool(), false, true, "BETA"},
                {new TestProductionTool(), false, false, null}
        };
    }

    @Test(dataProvider = "toolMaturityLevels")
    public void testToolMaturityLevel(
            final CommandLineProgram clp,
            final boolean isExperimental,
            final boolean isBeta,
            final String warningSentinel)
    {
        Assert.assertEquals(clp.isBetaFeature(), isBeta);
        Assert.assertEquals(clp.isExperimentalFeature(), isExperimental);
        if (warningSentinel == null) {
            Assert.assertEquals(clp.getToolStatusWarning(true), null);
        } else {
            Assert.assertTrue(clp.getToolStatusWarning(true).contains(warningSentinel));
        }
    }

}