package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.tribble.Tribble;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class VcfFormatConverterTest extends CommandLineProgramTest {
    private static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/vcf/vcfFormatTest");
    private static final String TEST_FILE_BASE = "vcfFormatTest";

    private static final String VCF = ".vcf";
    private static final String VCF_GZ = ".vcf.gz";
	private static final String BCF = ".bcf";

    private static final File TEST_VCF = new File(TEST_DATA_PATH, TEST_FILE_BASE + VCF);
    private static final File TEST_BCF = new File(TEST_DATA_PATH, TEST_FILE_BASE + BCF);

    public String getTestedClassName() {
        return VcfFormatConverter.class.getSimpleName();
    }

    @Test
    public void testVcfToVcf() {
        runLikeTest(TEST_VCF, VCF);
    }

    @Test
    public void testVcfToBcf() {
        runBackAndForthTest(TEST_VCF, BCF, VCF);
    }

    @Test
    public void testVcfToVcfGz() {
        runBackAndForthTest(TEST_VCF, VCF_GZ, VCF);
    }

    @Test
    public void testBcfToBcf() {
        runLikeTest(TEST_BCF, BCF);
    }

    @Test
    public void testBcfToVcf() {
        runBackAndForthTest(TEST_BCF, VCF, BCF);
    }

    private void runLikeTest(final File input, final String format) {
        final File outputFile = convertFile(input, "likeTest", format);
        compareFiles(input, outputFile);
    }

    private void runBackAndForthTest(final File input, final String format, final String originalFormat) {
        final String tempPrefix = "backAndForth";

        final File backAndForth = convertFile(input, tempPrefix, format);
        final File backAndForthSeries2 = convertFile(backAndForth, tempPrefix, originalFormat);

        compareFiles(input, backAndForthSeries2);
    }

    private File convertFile(final File input, final String prefix, final String format) {
        final File outputFile;
        try {
            outputFile = File.createTempFile(prefix, format);
        } catch (final IOException ioe) {
            throw new GATKException("Unable to create temp file!");
        }

        outputFile.deleteOnExit();
        new File(outputFile.getAbsolutePath() + Tribble.STANDARD_INDEX_EXTENSION).deleteOnExit();

        final List<String> args = new ArrayList<>();
        args.add("--INPUT");
        args.add(input.getAbsolutePath());
        args.add("--OUTPUT");
        args.add(outputFile.getAbsolutePath());
        if (VCF_GZ.equals(format)) {
            args.add("--CREATE_INDEX");
            args.add("false");
        }
        if (input.getName().endsWith(VCF_GZ)) {
            args.add("--REQUIRE_INDEX");
            args.add("false");
        }
        runCommandLine(args);
        return outputFile;
    }

    private void compareFiles(final File file1, final File file2) {
        // Ok, so this isn't exactly comparing md5 checksums or anything, but it should be good enough
        // for our purposes.
        Assert.assertTrue(file1.exists());
        Assert.assertTrue(file2.exists());
        Assert.assertEquals(file1.length(), file2.length());
    }

}
