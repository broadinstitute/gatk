package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.RevertSam;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class RevertSamSparkUnitTest extends CommandLineProgramTest {

    private final File basicSamToRevert = getTestFile("revert_sam_basic.sam");

    @DataProvider(name="positiveTestData")
    public Object[][] getPostitiveTestData() {
        return new Object[][] {
                {null, true, true, true, true, null, null, Collections.EMPTY_LIST},
                {SAMFileHeader.SortOrder.queryname, true, true, true, false, "Hey,Dad!", null, Arrays.asList("XT")},
                {null, false, true, false, false, "Hey,Dad!", "NewLibraryName", Arrays.asList("XT")},
                {null, false, false, false, false, null, null, Collections.EMPTY_LIST}
        };
    }

    @Test(dataProvider="positiveTestData")
    public void basicPositiveTests(final SAMFileHeader.SortOrder so, final boolean removeDuplicates, final boolean removeAlignmentInfo,
                                   final boolean restoreOriginalQualities, final boolean outputByReadGroup, final String sample, final String library,
                                   final List<String> attributesToClear) throws Exception {

        final File output = outputByReadGroup?Files.createTempDirectory("picardRevertSamTest").toFile():File.createTempFile("reverted", ".sam");
        File output0 = createTempFile("0", ".sam");
        File output1 = createTempFile("1", ".sam");
        File output2 = createTempFile("2", ".sam");
//        if (outputByReadGroup) {
//            output = Files.createTempDirectory("picardRevertSamTest").toFile();
//            output0 = Paths.get(output.toString(), "0.sam").toFile();
//            output1 = Paths.get(output.toString(), "1.sam").toFile();
//            output2 = Paths.get(output.toString(), "2.sam").toFile();
//        } else {
//            output = File.createTempFile("reverted", ".sam");
//        }
        output.deleteOnExit();
        final RevertSam reverter = new RevertSam();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(basicSamToRevert);
        args.addOutput(output);

        if (outputByReadGroup) {
            args.add("output-by-readgroup");
        }
        if (so != null) {
            args.addArgument("sort-order",so.name()); //TODO decide on sort order outputing
        }
//        args[index++] = "dontRemoveDuplicateInformation=" + removeDuplicates; //TODO this is unsuported
        args.add("remove-alignment-inormation");
        args.add("restore-original-qualities");
        if (sample != null) {
            args.addArgument("sample-alias",sample);
        }
        if (library != null) {
            args.addArgument("library-name",library);
        }
        for (final String attr : attributesToClear) {
            args.addArgument("attributes-to-clear",attr);
        }

        runCommandLine(args);

//        if (outputByReadGroup) {
//            verifyPositiveResults(output0, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "0", 2, sample, library);
//            verifyPositiveResults(output1, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "1", 4, sample, library);
//            verifyPositiveResults(output2, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "2", 2, sample, library);
//        } else {
//            verifyPositiveResults(output, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, null, 8, sample, library);
//        }
    }
}
