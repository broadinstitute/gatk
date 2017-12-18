package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.PrintReads;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Tests for validation of sequence dictionaries.
 */
public class SequenceDictionaryValidationWalkerIntegrationTest extends CommandLineProgramTest {

    public static final String SEQDICTVAL_TEST_DIRECTORY = "src/test/resources/org/broadinstitute/hellbender/utils/SequenceDictionaryUtils/";

    @Override
    public String getTestedClassName() {
        return PrintReads.class.getSimpleName();
    }

    @Test
    public void testSeqDictValCompatible() throws IOException {
        final File outFile = GATKBaseTest.createTempFile("testSeqDictValCompatible", ".bam");
        final String[] args = {
                "--input" , SEQDICTVAL_TEST_DIRECTORY + "test.sorted.bam",
                "-R", SEQDICTVAL_TEST_DIRECTORY + "test.fasta",
                "--output", outFile.getAbsolutePath()
        };
        // Should run without producing an exception, since the sequence dictionaries are compatible
        runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testSeqDictValIncompatible() throws IOException {
        final File outFile = GATKBaseTest.createTempFile("testSeqDictValIncompatible", ".bam");
        final String[] args = {
                "--input" , SEQDICTVAL_TEST_DIRECTORY + "test2.sorted.bam",
                "-R", SEQDICTVAL_TEST_DIRECTORY + "test2.fasta",
                "--output", outFile.getAbsolutePath()
        };
        // Should produce an exception, since the sequence dictionaries are incompatible
        runCommandLine(args);
    }

    @Test
    public void testSeqDictValIncompatibleDisableValidation() throws IOException {
        final File outFile = GATKBaseTest.createTempFile("testSeqDictValIncompatibleDisableValidation", ".bam");
        final String[] args = {
                "--input" , SEQDICTVAL_TEST_DIRECTORY + "test2.sorted.bam",
                "-R", SEQDICTVAL_TEST_DIRECTORY + "test2.fasta",
                "--output", outFile.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME
        };
        // Should run without producing an exception even though the dictionaries are incompatible,
        // since we've disabled sequence dictionary validation
        runCommandLine(args);
    }
}
