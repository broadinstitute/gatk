package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public final class MergeSamFilesIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/MergeSamFiles");

    public String getTestedClassName() {
        return MergeSamFiles.class.getSimpleName();
    }
    /**
     * Confirm that unsorted input can result in coordinate sorted output, with index created.
     */
    @Test
    public void unsortedInputSortedOutputTest() throws Exception {
        final File unsortedInputTestDataDir = new File(TEST_DATA_DIR, "unsorted_input");
        final File sam1 = new File(unsortedInputTestDataDir, "1.sam");
        final File sam2 = new File(unsortedInputTestDataDir, "2.sam");
        final File mergedOutput = BaseTest.createTempFile("unsortedInputSortedOutputTest.", BamFileIoUtils.BAM_FILE_EXTENSION);

        final String[] args = new String[]{
                "--input", sam1.getAbsolutePath(),
                "--input", sam2.getAbsolutePath(),
                "--output", mergedOutput.getAbsolutePath(),
                "--SO", "coordinate"
        };
        
        runCommandLine(args);
        final SamReader reader = SamReaderFactory.makeDefault().open(mergedOutput);
        Assert.assertEquals(reader.getFileHeader().getSortOrder(), SAMFileHeader.SortOrder.coordinate);
        SamAssertionUtils.assertSamValid(mergedOutput);
        CloserUtil.close(reader);
    }
}
