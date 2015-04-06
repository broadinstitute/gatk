package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author George Grant
 */
public class UpdateVcfSequenceDictionaryTest extends CommandLineProgramTest {
    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/vcf/vcfFormatTest");
    private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("UpdateVcfSequenceDictionaryTest", null);

    @AfterClass
    public void teardown() {
        IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
    }

    @Test
    public void testUpdateVcfSequenceDictionary() {
        final File input = new File(TEST_DATA_PATH, "vcfFormatTest.vcf");
        // vcfFormatTest.bad_dict.vcf is a vcf with two (2) ##contig lines deleted
        final File samSequenceDictionaryVcf = new File(TEST_DATA_PATH, "vcfFormatTest.bad_dict.vcf");
        final File outputFile = new File(OUTPUT_DATA_PATH, "updateVcfSequenceDictionaryTest-delete-me.vcf");

        outputFile.deleteOnExit();

        final String[] argv = new String[]{
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outputFile.getAbsolutePath(),
                "--SEQUENCE_DICTIONARY", samSequenceDictionaryVcf.getAbsolutePath()
        };

        runCommandLine(argv);

        IOUtil.assertFilesEqual(samSequenceDictionaryVcf, outputFile);

        // A little extra checking.
        Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(input).size(), 3);
        Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(samSequenceDictionaryVcf).size(), 1);
        Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(outputFile).size(), 1);
    }
}
