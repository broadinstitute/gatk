package org.broadinstitute.hellbender.testutils;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class SamAssertionUtilsUnitTest extends GATKBaseTest {

    @DataProvider(name = "bamPairs")
    public Object[][] differentFilesButSameContent(){
        final String fileBam= "file1.bam";
        final String fileSam = "file1.sam";
        return new Object[][]{
                {fileBam, fileBam, true},
                {fileSam, fileSam, true},
                {fileBam, fileSam, true},
                {fileSam, fileBam, true},
                {fileBam, "file1_reorder_PG_header_lines.sam", true},    //ok to reorder PG lines
                {fileBam, "file1_reorder_read_attributes.sam", true},    //ok to reorder attributes
                {fileBam, "file1_new_read_attribute.sam", true},         //ok to add an attribute
                {fileBam, "file1_different_version.sam", true},         //ok to have a different version if everything else is same

                {fileBam, "file1_missing_read_attribute.sam", false},  //not ok to lose an attribute
                {fileBam, "file1_different_bases.sam", false},
                {fileBam, "file1_different_basequals.sam", false},
                {fileBam, "file1_different_attributes.sam", false},
                {fileBam, "file1_different_mappingQ.sam", false},
                {fileBam, "file1_different_cigar.sam", false},
                {fileBam, "file1_different_position.sam", false},
        };
    }

    private static final File TEST_DATA_DIR = new File(publicTestDir, "org/broadinstitute/hellbender/utils/test/SamAssertionUtilsUnitTest");

    @Test(dataProvider = "bamPairs")
    public void testCompareStringent(final String expectedFile, final String actualFile, boolean expectedEqual) throws Exception {
        final File expectedF= new File(TEST_DATA_DIR, expectedFile);
        final File actualF= new File(TEST_DATA_DIR, actualFile);
        Assert.assertEquals(expectedEqual, null == SamAssertionUtils.samsEqualStringent(actualF, expectedF, ValidationStringency.LENIENT, null));
    }

    @DataProvider(name="testCRAMContentsSucceed")
    public Object[][] testCRAMContentsSucceedData() {
        return new Object[][] {
                {new File(TEST_DATA_DIR, "valid.cram")},
        };
    }

    @DataProvider(name="testCRAMContentsFail")
    public Object[][] testCRAMContentsFailData() {
        return new Object[][] {
                {new File(TEST_DATA_DIR, "file1.sam") },
                {new File(TEST_DATA_DIR, "file1.bam") },
                {new File(TEST_DATA_DIR, "fake_cram_with_bam_contents.cram") }
        };
    }

    @Test(dataProvider = "testCRAMContentsSucceed")
    public void testAssertCRAMContentsSucceed(File putativeCRAMFile) {
        SamAssertionUtils.assertCRAMContents(putativeCRAMFile.toPath());
    }

    @Test(dataProvider = "testCRAMContentsFail", expectedExceptions=AssertionError.class)
    public void testAssertCRAMContentsFail(File putativeCRAMFile) {
        SamAssertionUtils.assertCRAMContents(putativeCRAMFile.toPath());
    }

    @DataProvider(name="testCRAMContentsIfCRAMSucceed")
    public Object[][] testCRAMContentsIfCRAMSucceedData() {
        return new Object[][] {
                {new File(TEST_DATA_DIR, "file1.sam") },
                {new File(TEST_DATA_DIR, "file1.bam") },
                {new File(TEST_DATA_DIR, "valid.cram") },
        };
    }

    @DataProvider(name="testCRAMContentsIfCRAMFail")
    public Object[][] testCRAMContentsIfCRAMFailData() {
        return new Object[][] {
                {new File(TEST_DATA_DIR, "fake_cram_with_bam_contents.cram") }
        };
    }

    @Test(dataProvider = "testCRAMContentsIfCRAMSucceed")
    public void testAssertCRAMContentsIfCRAMSucceed(File putativeCRAMFile) {
        SamAssertionUtils.assertCRAMContentsIfCRAM(putativeCRAMFile);
    }

    @Test(dataProvider = "testCRAMContentsIfCRAMFail", expectedExceptions = AssertionError.class)
    public void testAssertCRAMContentsIfCRAM(File putativeCRAMFile) {
        SamAssertionUtils.assertCRAMContentsIfCRAM(putativeCRAMFile);
    }
}
