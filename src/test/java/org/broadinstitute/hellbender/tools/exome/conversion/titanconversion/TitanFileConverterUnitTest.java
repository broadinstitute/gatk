package org.broadinstitute.hellbender.tools.exome.conversion.titanconversion;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;


public class TitanFileConverterUnitTest extends BaseTest {
    private static final String ACNV_TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome";

    // Typically, this next file would be the output of a tangent normalization run.
    private static final File COVERAGES_FILE = new File(ACNV_TEST_SUB_DIR, "coverages-for-allelic-integration.tsv");
    private static final File TUMOR_ALLELIC_COUNTS_FILE = new File(ACNV_TEST_SUB_DIR, "snps-full.tsv");

    @Test
    public void basicTNConverterTest() {
        final File tempOutput = createTempFile("titanTN", ".tsv");
        TitanFileConverter.convertCRToTitanCovFile(COVERAGES_FILE, tempOutput);
        try {
            Assert.assertEquals(FileUtils.readLines(tempOutput).size(), FileUtils.readLines(COVERAGES_FILE).size());
            final String[] headers = ArrayUtils.toArray(FileUtils.readLines(tempOutput).get(0).split("\t"));
            Assert.assertEquals(headers, TitanCopyRatioEstimateColumns.FULL_COLUMN_NAME_ARRAY.toArray());
        } catch (final IOException ioe) {
            Assert.fail("Problem with unit test configuration.", ioe);
        }
    }
    @Test
    public void basicHetConverterTest() {
        final File tempOutput = createTempFile("titanHet", ".tsv");
        TitanFileConverter.convertHetPulldownToTitanHetFile(TUMOR_ALLELIC_COUNTS_FILE, tempOutput);
        try {
            Assert.assertEquals(FileUtils.readLines(tempOutput).size(), FileUtils.readLines(TUMOR_ALLELIC_COUNTS_FILE).size());
            final List<String> headers = Arrays.asList(FileUtils.readLines(tempOutput).get(0).split("\t"));
            Assert.assertEquals(headers, TitanAllelicCountTableColumn.COLUMNS.names());
        } catch (final IOException ioe) {
            Assert.fail("Problem with unit test configuration.", ioe);
        }
    }
}
