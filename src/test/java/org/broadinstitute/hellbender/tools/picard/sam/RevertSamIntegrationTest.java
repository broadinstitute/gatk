package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public final class RevertSamIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = getTestDataDir();
    public static final File basicSamToRevert = new File(TEST_DATA_DIR, "revert_sam_basic.sam");
    public static final File negativeTestSamToRevert = new File(TEST_DATA_DIR, "revert_sam_negative.sam");

    private static final String revertedQualities  =
            "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111";

    private static final String unmappedRead = "both_reads_present_only_first_aligns/2";

    @Override
    public String getTestedClassName() {
        return RevertSam.class.getSimpleName();
    }

    @Test(dataProvider="positiveTestData")
    public void basicPositiveTests(final SAMFileHeader.SortOrder so, final boolean removeDuplicates, final boolean removeAlignmentInfo,
                                   final boolean restoreOriginalQualities, final String sample, final String library,
                                   final List<String> attributesToClear) throws Exception {

        final File output = BaseTest.createTempFile("reverted", ".sam");
        final RevertSam reverter = new RevertSam();
        final List<String> args = new ArrayList<>();
        args.add("--input");
        args.add(basicSamToRevert.getAbsolutePath());
        args.add("--output");
        args.add(output.getAbsolutePath());
        if (so != null) {
            args.add("--SORT_ORDER");
            args.add(so.name());
        }
        args.add("--REMOVE_DUPLICATE_INFORMATION");
        args.add(Boolean.toString(removeDuplicates));
        args.add("--REMOVE_ALIGNMENT_INFORMATION");
        args.add(Boolean.toString(removeAlignmentInfo));
        args.add("--RESTORE_ORIGINAL_QUALITIES");
        args.add(Boolean.toString(restoreOriginalQualities));
        if (sample != null) {
            args.add("--SAMPLE_ALIAS");
            args.add(sample);
        }
        if (library != null) {
            args.add("--LIBRARY_NAME");
            args.add(library);
        }
        for (final String attr : attributesToClear) {
            args.add("--ATTRIBUTE_TO_CLEAR");
            args.add(attr);
        }
        runCommandLine(args);

        final SamReader reader = SamReaderFactory.makeDefault().open(output);
        final SAMFileHeader header = reader.getFileHeader();
        Assert.assertEquals(header.getSortOrder(), SAMFileHeader.SortOrder.queryname);
        Assert.assertEquals(header.getProgramRecords().size(), removeAlignmentInfo ? 0 : 1);
        for (final SAMReadGroupRecord rg : header.getReadGroups()) {
            if (sample != null) {
                Assert.assertEquals(rg.getSample(), sample);
            }
            else {
                Assert.assertEquals(rg.getSample(), "Hi,Mom!");
            }
            if (library != null) {
                Assert.assertEquals(rg.getLibrary(), library);
            }
            else {
                Assert.assertEquals(rg.getLibrary(), "my-library");
            }
        }
        for (final SAMRecord rec : reader) {
            if (removeDuplicates) {
                Assert.assertFalse(rec.getDuplicateReadFlag(),
                        "Duplicates should have been removed: " + rec.getReadName());
            }

            if (removeAlignmentInfo) {
                Assert.assertTrue(rec.getReadUnmappedFlag(),
                        "Alignment info should have been removed: " + rec.getReadName());
            }

            if (restoreOriginalQualities && !unmappedRead.equals(
                    rec.getReadName() + "/" + (rec.getFirstOfPairFlag() ? "1" : "2"))) {

                Assert.assertEquals(rec.getBaseQualityString(), revertedQualities);
            } else {
                Assert.assertNotSame(rec.getBaseQualityString(), revertedQualities);
            }

            for (final SAMRecord.SAMTagAndValue attr : rec.getAttributes()) {
                if (removeAlignmentInfo || (!attr.tag.equals("PG") && !attr.tag.equals("NM")
                        && !attr.tag.equals("MQ"))) {
                    Assert.assertFalse(reverter.ATTRIBUTE_TO_CLEAR.contains(attr.tag),
                            attr.tag + " should have been cleared.");
                }
            }
        }
        CloserUtil.close(reader);
    }


    @DataProvider(name="positiveTestData")
    public Object[][] getPostitiveTestData() {
        return new Object[][] {
                {null, true, true, true, null, null, Collections.EMPTY_LIST},
                {SAMFileHeader.SortOrder.queryname, true, true, true, "Hey,Dad!", null,
                        Arrays.asList("XT", SAMTag.NM.name(), SAMTag.UQ.name(), SAMTag.PG.name(), SAMTag.MD.name(),
                                SAMTag.MQ.name(), SAMTag.SA.name(), SAMTag.MC.name())},
                {null, false, true, false, "Hey,Dad!", "NewLibraryName",
                        Arrays.asList("XT", SAMTag.NM.name(), SAMTag.UQ.name(), SAMTag.PG.name(), SAMTag.MD.name(),
                        SAMTag.MQ.name(), SAMTag.SA.name(), SAMTag.MC.name())},
                {null, false, false, false, null, null, Collections.EMPTY_LIST}
        };
    }


    @Test(dataProvider="negativeTestData", expectedExceptions = {UserException.class, GATKException.class})
    public void basicNegativeTest(final String sample, final String library) throws Exception {

        final File output = BaseTest.createTempFile("bad", ".sam");
        final RevertSam reverter = new RevertSam();
        final List<String> args = new ArrayList<>();
        args.add("--input");
        args.add(negativeTestSamToRevert.getAbsolutePath());
        args.add("--output");
        args.add(output.getAbsolutePath());
        if (sample != null) {
            args.add("--SAMPLE_ALIAS");
            args.add(sample);
        }
        if (library != null) {
            args.add("--LIBRARY_NAME");
            args.add(library);
        }
        runCommandLine(args);
        Assert.fail("Negative test should have thrown an exception and didn't");
    }

    @DataProvider(name="negativeTestData")
    public Object[][] getNegativeTestData() {
        return new Object[][] {
                {"NewSample", null},
                {null, "NewLibrary"},
                {"NewSample", "NewLibrary"}
        };
    }
}
