package org.broadinstitute.hellbender.tools.picard;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

public class MergeSamFilesTest extends CommandLineProgramTest {
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
        assertSamValid(mergedOutput);
        CloserUtil.close(reader);
    }

    private void assertSamValid(final File sam) throws IOException {
        final SamReader samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(sam);
        final SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 8000);
        validator.setIgnoreWarnings(true);
        validator.setVerbose(true, 1000);
        validator.setErrorsToIgnore(Arrays.asList(SAMValidationError.Type.MISSING_READ_GROUP));
        final boolean validated = validator.validateSamFileVerbose(samReader, null);
        samReader.close();
        Assert.assertTrue(validated, "ValidateSamFile failed");
    }
}
