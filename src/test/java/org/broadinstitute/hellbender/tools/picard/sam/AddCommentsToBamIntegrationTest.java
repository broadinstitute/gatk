package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public final class AddCommentsToBamIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = getTestDataDir();
    private static final File BAM_FILE = new File(TEST_DATA_DIR, "add_comments_to_bam.bam");
    private static final File SAM_FILE = new File(TEST_DATA_DIR, "add_comments_to_bam.sam");

    private static final String[] commentList = new String[]{"test1", "test2", "test3"};

    public String getTestedClassName() {
        return AddCommentsToBam.class.getSimpleName();
    }

    @Test
    public void testAddCommentsToBam() throws Exception {
        final File outputFile = BaseTest.createTempFile("addCommentsToBamTest.", BamFileIoUtils.BAM_FILE_EXTENSION);
        runIt(BAM_FILE, outputFile, commentList);

        final SAMFileHeader newHeader = SamReaderFactory.makeDefault().getFileHeader(outputFile);

        // The original comments are massaged when they're added to the header. Perform the same massaging here,
        // and then compare the lists
        final List<String> massagedComments = new LinkedList<>();
        for (final String comment : commentList) {
            massagedComments.add(SAMTextHeaderCodec.COMMENT_PREFIX + comment);
        }

        Assert.assertEquals(newHeader.getComments(), massagedComments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testUsingSam() throws Exception {
        final File outputFile = BaseTest.createTempFile("addCommentsToBamTest.samFile", BamFileIoUtils.BAM_FILE_EXTENSION);
        runIt(SAM_FILE, outputFile, commentList);
        throw new IllegalStateException("We shouldn't be here!");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testUsingNewlines() throws Exception {
        final File outputFile = BaseTest.createTempFile("addCommentsToBamTest.newLine", BamFileIoUtils.BAM_FILE_EXTENSION);
        runIt(SAM_FILE, outputFile, new String[]{"this is\n a crazy\n test"});
        throw new IllegalStateException("We shouldn't be here!");
    }

    private void runIt(final File inputFile, final File outputFile, final String[] commentList) {
        final List<String> args = new ArrayList<>(Arrays.asList(
                "--input", inputFile.getAbsolutePath(),
                "--output", outputFile.getAbsolutePath()));
        for (final String comment : commentList) {
            args.add("--COMMENT");
            args.add(comment);
        }
        runCommandLine(args);
    }

}
