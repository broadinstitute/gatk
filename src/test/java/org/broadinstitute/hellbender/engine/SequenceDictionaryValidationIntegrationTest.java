package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;

/**
 * Created by edwardk on 7/6/15.
 */
public class SequenceDictionaryValidationIntegrationTest extends CommandLineProgramTest {

    public static final String SEQDICTVAL_TEST_DIRECTORY = "src/test/resources/org/broadinstitute/hellbender/utils/SequenceDictionaryUtils/";

    @Override
    public String getTestedClassName() {
        return "ExampleIntervalWalker";
    }

    @Test
    public void testSeqDictValCompatible() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + SEQDICTVAL_TEST_DIRECTORY + "test.fasta" +
                        " -I " + SEQDICTVAL_TEST_DIRECTORY + "test.sorted.bam" +
                        " -V " + "TestFeatures:" + SEQDICTVAL_TEST_DIRECTORY + "test.vcf" +
                        " -L " + SEQDICTVAL_TEST_DIRECTORY + "test.intervals",
                Collections.emptyList());
        testSpec.executeTest("testSequenceDictionariesCompatible", this);
    }

    @Test
    public void testSeqDictValIncompatible() throws IOException {
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -R " + SEQDICTVAL_TEST_DIRECTORY + "test2.fasta" +
                        " -I " + SEQDICTVAL_TEST_DIRECTORY + "test2.sorted.bam" +
                        " -V " + "TestFeatures:" + SEQDICTVAL_TEST_DIRECTORY + "test2.vcf" +
                        " -L " + SEQDICTVAL_TEST_DIRECTORY + "test2.intervals",
                0,
                UserException.IncompatibleSequenceDictionaries.class);
        testSpec.executeTest("testSequenceDictionariesIncompatible", this);
    }
}
