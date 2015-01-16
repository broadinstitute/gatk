package org.broadinstitute.hellbender.tools.picard;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.sam.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class MergeSamFilesIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/MergeSamFiles");

    public String getCommandLineProgramName() {
        return MergeSamFiles.class.getSimpleName();
    }
    /**
     * Confirm that unsorted input can result in coordinate sorted output, with index created.
     */
    @Test
    public void unsortedInputSortedOutputTest() throws Exception {
        final File unsortedInputTestDataDir = new File(TEST_DATA_DIR, "unsorted_input");
        final File mergedOutput = File.createTempFile("unsortedInputSortedOutputTest.", BamFileIoUtils.BAM_FILE_EXTENSION);
        mergedOutput.deleteOnExit();
        final ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("I=" + new File(unsortedInputTestDataDir, "1.sam").getAbsolutePath());
        args.add("I=" + new File(unsortedInputTestDataDir, "2.sam").getAbsolutePath());
        args.add("O=" + mergedOutput.getAbsolutePath());
        args.add("SO=coordinate");

        Assert.assertEquals(runCommandLine(args.getArgsList()), null);
        final SamReader reader = SamReaderFactory.makeDefault().open(mergedOutput);
        Assert.assertEquals(reader.getFileHeader().getSortOrder(), SAMFileHeader.SortOrder.coordinate);
        SamAssertionUtils.assertSamValid(mergedOutput);
        CloserUtil.close(reader);
    }
}
