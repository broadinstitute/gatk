package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.RevertSam;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

@Test(groups = "Spark")
public class RevertSamSparkIntegrationTest extends CommandLineProgramTest {

    private static final List<String> defaultAttributesToClearPlusXT = new ArrayList<String>() {
        private static final long serialVersionUID = 1L;
        {
        addAll(RevertSamSpark.DEFAULT_ATTRIBUTES_TO_CLEAR);
        add("XT");
    }};

    private final File basicSamToRevert = getTestFile("revert_sam_basic.sam");
    private final File basicBamToRevert = getTestFile("revert_sam_basic.bam");
    private final File basicCramToRevert = getTestFile("revert_sam_basic.cram");
    private final File sampleLibraryOverrideSam = getTestFile("revert_sam_sample_library_override.sam");
    private final File validOutputMap = getTestFile("revert_sam_valid_output_map.txt");
    private final File nonExistentOutputMap = getTestFile("revert_sam_does_not_exist.txt");
    private final File badHeaderOutputMap = getTestFile("revert_sam_bad_header_output_map.txt");
    private final File referenceFasta = getTestFile("test.fasta");
    private final File singleEndSamToRevert = getTestFile("revert_sam_single_end.sam");
    private final File missingRGInfo = getTestFile("missing-rg-info.sam");

    private static final String revertedQualities  =
            "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111";

    private static final String unmappedRead = "both_reads_present_only_first_aligns/2";

    @DataProvider(name="positiveTestData")
    public Object[][] positiveTestData() {
        return new Object[][] {
                {basicSamToRevert, null, false, false, true, true, null, null, Collections.EMPTY_LIST},
                {basicSamToRevert, SAMFileHeader.SortOrder.queryname, false, false, true, false, "Hey,Dad!", null, defaultAttributesToClearPlusXT},
                {basicSamToRevert, null, true, false, false, false, "Hey,Dad!", "NewLibraryName", defaultAttributesToClearPlusXT},
                {basicSamToRevert, null, true, true, false, false, null, null, Collections.EMPTY_LIST},

                {basicBamToRevert, null, false, false, true, true, null, null, Collections.EMPTY_LIST},
                {basicBamToRevert, SAMFileHeader.SortOrder.queryname, false, false, true, false, "Hey,Dad!", null, defaultAttributesToClearPlusXT},
                {basicBamToRevert, null, true, false, false, false, "Hey,Dad!", "NewLibraryName", defaultAttributesToClearPlusXT},
                {basicBamToRevert, null, true, true, false, false, null, null, Collections.EMPTY_LIST},

                {basicCramToRevert, null, false, false, true, true, null, null, Collections.EMPTY_LIST},
                {basicCramToRevert, SAMFileHeader.SortOrder.queryname, false, false, true, false, "Hey,Dad!", null, defaultAttributesToClearPlusXT},
                {basicCramToRevert, null, true, false, false, false, "Hey,Dad!", "NewLibraryName", defaultAttributesToClearPlusXT},
                {basicCramToRevert, null, true, true, false, false, null, null, Collections.EMPTY_LIST}
        };
    }

    @Test(dataProvider= "positiveTestData")
    public void basicPositiveTests(final File input, final SAMFileHeader.SortOrder so, final boolean removeDuplicates, final boolean removeAlignmentInfo,
                                   final boolean restoreOriginalQualities, final boolean outputByReadGroup, final String sample, final String library,
                                   final List<String> attributesToClear) throws Exception {

        final File output = outputByReadGroup ? Files.createTempDirectory("picardRevertSamSparkTest").toFile() : BaseTest.createTempFile("reverted", ".sam");
        final File output0 = new File(output.getPath()+"/0.sam");
        final File output1 = new File(output.getPath()+"/1.sam");
        final File output2 = new File(output.getPath()+"/2.sam");
        final File output3 = new File(output.getPath()+"/3.sam");
        final RevertSamSpark reverter = new RevertSamSpark();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(input);
        args.addOutput(output);
        args.addReference(referenceFasta);

        if (outputByReadGroup) {
            args.add(RevertSamSpark.OUTPUT_BY_READGROUP_FILE_FORMAT_LONG_NAME, RevertSamSpark.FileType.sam);
            args.add(RevertSamSpark.OUTPUT_BY_READGROUP_LONG_NAME, true);
        }
        if (so != null) {
            args.add("sort-order", so); //TODO decide on sort order outputing
        }
        if (!removeAlignmentInfo) {
            args.add(RevertSamSpark.KEEP_ALIGNMENT_INFORMATION, true);
        }
        if (sample != null) {
            args.add("sample-alias", sample);
        }
        if (library != null) {
            args.add("library-name", library);
        }
        for (final String attr : attributesToClear) {
            args.add("attributes-to-clear", attr);
        }

        runCommandLine(args);

        if (outputByReadGroup) {
            verifyPositiveResults(output0, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "0", 2, sample, library);
            verifyPositiveResults(output1, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "1", 4, sample, library);
            verifyPositiveResults(output2, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "2", 2, sample, library);
            verifyPositiveResults(output3, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, "3", 0, sample, library);
        } else {
            verifyPositiveResults(output, reverter, removeDuplicates, removeAlignmentInfo, restoreOriginalQualities, outputByReadGroup, null, 8, sample, library);
        }
    }

    @DataProvider
    public Object[][] getInputs(){
        return new Object[][] {
                {basicSamToRevert},
                {basicBamToRevert},
                {basicCramToRevert}
        };
    }

    @Test(dataProvider = "getInputs")
    public void testOutputByReadGroupWithOutputMap(File input) throws Exception {
        final File outputDir = createTempDir("testOutputByReadGroupWithOutputMap");
        // Create the output map
        final File outputMapFile = BaseTest.createTempFile("picardRevertSamSparkTestOutputMap", ".txt");
        final String outputPath0;
        final String outputPath1;
        final String outputPath2;
        try (final PrintWriter mapWriter = new PrintWriter(outputMapFile)) {
            outputPath0 = outputDir + "/my_rg0.sam";
            outputPath1 = outputDir + "/rg1.sam";
            outputPath2 = outputDir + "/my_rg2.bam";
            final String outputPath3 = outputDir + "/my_rg3.sam";//TODO not used?
            mapWriter.println("READ_GROUP_ID\tOUTPUT");
            mapWriter.println("0\t" + outputPath0);
            mapWriter.println("2\t" + outputPath2);
            mapWriter.println("1\t" + outputPath1);
            mapWriter.println("3\t" + outputPath3);
            Assert.assertFalse(mapWriter.checkError(), "PrintWriter encountered an error");
        }

        final RevertSamSpark reverter = new RevertSamSpark();

        final String args[] = new String[] {
                "-I",input.getPath(),
                "--output-by-readgroup",
                "--output-by-readgroup-file-format", "sam",
                "--output-map",outputMapFile.getPath(),
                "-R",referenceFasta.getPath(),
                "--sort-order",SAMFileHeader.SortOrder.queryname.name(),
                "--"+RevertSamSpark.SAMPLE_ALIAS_ARG,"test_sample_1",
                "--"+RevertSamSpark.LIBRARY_NAME_ARG,"test_library_1",
                "--"+RevertSamSpark.ATTRIBUTE_TO_CLEAR_LONG_NAME,SAMTag.NM.name()
        };

        runCommandLine(args);

        final File output0 = new File(outputPath0);
        final File output1 = new File(outputPath1);
        final File output2 = new File(outputPath2);
        verifyPositiveResults(output0, reverter, true, true, true, true, "0", 2, "test_sample_1", "test_library_1");
        verifyPositiveResults(output1, reverter, true, true, true, true, "1", 4, "test_sample_1", "test_library_1");
        verifyPositiveResults(output2, reverter, true, true, true, true, "2", 2, "test_sample_1", "test_library_1");
    }

    @Test(dataProvider = "getInputs")
    public void testSingleEndSanitize(File input) throws Exception {
        final File output = createTempFile("single_end_reverted", ".sam");
        final String[] args = { "-I", input.getAbsolutePath(), "-O", output.getAbsolutePath(), "--sanitize", "-R", referenceFasta.getAbsolutePath()};
        runCommandLine(args);
    }

    private void verifyPositiveResults(
            final File outputFile,
            final RevertSamSpark reverter,
            final boolean removeDuplicates,
            final boolean removeAlignmentInfo,
            final boolean restoreOriginalQualities,
            final boolean outputByReadGroup,
            final String readGroupId,
            final int numReadsExpected,
            final String sample,
            final String library) throws IOException {

        try (SamReader reader = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(outputFile)) {
            final SAMFileHeader header = reader.getFileHeader();
            Assert.assertEquals(header.getSortOrder(), SAMFileHeader.SortOrder.queryname);
            Assert.assertEquals(header.getProgramRecords().size(), removeAlignmentInfo ? 0 : 1);
            final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
            if (outputByReadGroup) {
                Assert.assertEquals(readGroups.size(), 1);
                Assert.assertEquals(readGroups.get(0).getId(), readGroupId);
            }
            for (final SAMReadGroupRecord rg : header.getReadGroups()) {
                if (sample != null) {
                    Assert.assertEquals(rg.getSample(), sample);
                } else {
                    Assert.assertEquals(rg.getSample(), "Hi,Mom!");
                }
                if (library != null) {
                    Assert.assertEquals(rg.getLibrary(), library);
                } else {
                    Assert.assertEquals(rg.getLibrary(), "my-library");
                }
            }
            int numReads = 0;
            for (final SAMRecord rec : reader) {
                numReads++;
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
                            && !attr.tag.equals(SAMTag.MQ.toString()))) {
                        Assert.assertFalse(reverter.attributesToClear.contains(attr.tag),
                                attr.tag + " should have been cleared.");
                    }
                }
            }
            Assert.assertEquals(numReads, numReadsExpected);
        }
    }

    @DataProvider
    public Object[][] getFileTypes(){
        return new Object[][] {
                {".sam"},
                {".bam"},
                {".cram"}
        };
    }

    @Test(dataProvider = "getFileTypes")
    public void testSanitizeAndDeduplicateRecords(String extension) throws Exception {
        final File input  = BaseTest.createTempFile("test-input-santize-and-deduplicate-records", extension);
        final File output = BaseTest.createTempFile("test-output-santize-and-deduplicate-records", extension);

        // Create a SAM file that has duplicate records
        int numDuplicated;
        try (final SamReader reader = SamReaderFactory.makeDefault().open(basicSamToRevert);
             final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(reader.getFileHeader(), false, input, referenceFasta)) {

            numDuplicated = 0;
            for (final SAMRecord rec : reader) {
                writer.addAlignment(rec);
                if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) {
                    writer.addAlignment(rec);
                    numDuplicated++;
                }
            }
        }

        // Make sure some records are duplicated
        Assert.assertTrue(numDuplicated > 0);

        final String [] args = new String[]{
                "--input", input.getAbsolutePath(),
                "--sanitize",
                "--keep-first-duplicate",
                "--"+RevertSamSpark.DONT_RESTORE_ORIGINAL_QUALITIES_LONG_NAME,
                "-O", output.getAbsolutePath(),
                "-R", referenceFasta.getAbsolutePath()
        };
        runCommandLine(args);
        verifyPositiveResults(output, new RevertSamSpark(), false, true, false, false, null, 8, null, null);
    }

    @Test(dataProvider="overrideTestData", expectedExceptions = {UserException.class})
    public void testSampleLibraryOverride(final String sample, final String library) throws Exception {
        final File output = createTempFile("bad", ".sam");
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(sampleLibraryOverrideSam);
        args.addOutput(output);
        if (sample != null) {
            args.add(RevertSamSpark.SAMPLE_ALIAS_ARG,sample);
        }
        if (library != null) {
            args.add(RevertSamSpark.LIBRARY_NAME_ARG,library);
        }
        runCommandLine(args);
    }

    @DataProvider(name="overrideTestData")
    public Object[][] getNegativeTestData() {
        return new Object[][] {
                {"NewSample", null},
                {null, "NewLibrary"},
                {"NewSample", "NewLibrary"}
        };
    }

    @Test
    public void testIsOutputMapHeaderValid() {
        boolean isValid = RevertSamSpark.isOutputMapHeaderValid(Arrays.asList("READ_GROUP_ID","OUTPUT"));
        Assert.assertEquals(isValid, true);

        isValid = RevertSamSpark.isOutputMapHeaderValid(Arrays.asList("OUTPUT"));
        Assert.assertEquals(isValid, false);

        isValid = RevertSamSpark.isOutputMapHeaderValid(Collections.emptyList());
        Assert.assertEquals(isValid, false);
    }

    @Test
    public void testFilePathsWithMapFile() {
        final Map<String, Path> outputMap = RevertSamSpark.getOutputMap(validOutputMap.getAbsolutePath(), null, ".bam", Collections.emptyList(), true);
        Assert.assertEquals(outputMap.get("rg1"), IOUtils.getPath(new File("/path/to/my_rg_1.ubam").getAbsolutePath()));
        Assert.assertEquals(outputMap.get("rg2"), IOUtils.getPath(new File("/path/to/my_rg_2.ubam").getAbsolutePath()));
    }

    @Test
    public void testNoRgInfoSanitize() throws Exception {
        final File output = BaseTest.createTempFile("no-rg-reverted", ".sam");
        final String [] args = new String[]{
                "-I",missingRGInfo.getAbsolutePath(),
                "--sanitize",
                "-O", output.getAbsolutePath()
        };
        runCommandLine(args);
        verifyPositiveResults(output, new RevertSamSpark(), true,  true, false, false, null, 240, null, null);
    }

}