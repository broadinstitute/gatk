package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.SBIIndex;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;


public class CreateHadoopBamSplittingIndexIntegrationTest extends CommandLineProgramTest{

    private final File UNSORTED_BAM = getTestFile("count_reads.bam");
    private final File SORTED_BAM = getTestFile("count_reads_sorted.bam");

    @DataProvider(name="indexable")
    public Object[][] getIndexableFiles(){
        return new Object[][]{
                {UNSORTED_BAM},
                {SORTED_BAM}
        };
    }

    @Test(dataProvider = "indexable")
    public void testCreateSplittingIndex(final File bam) throws IOException {
        final File splittingIndex = getTempIndexFile();
        splittingIndex.delete();
        Assert.assertFalse(splittingIndex.exists());
        final ArgumentsBuilder args = getInputAndOutputArgs(bam, splittingIndex)
           .add("--"+ CreateHadoopBamSplittingIndex.SPLITTING_INDEX_GRANULARITY_LONG_NAME).add("1");
        this.runCommandLine(args);
        assertIndexIsNotEmpty(splittingIndex);

        //checked in index created with
        // ./gatk CreateHadoopBamSplittingIndex --input <filename> --splitting-index-granularity 1
        final File expectedSplittingIndex = new File(bam.toPath() + SBIIndex.FILE_EXTENSION);

        IOUtil.assertFilesEqual(splittingIndex, expectedSplittingIndex);
    }

    private static void assertIndexIsNotEmpty(final File splittingIndex) throws IOException {
        Assert.assertTrue(splittingIndex.exists());
        final SBIIndex splittingBAMIndex = SBIIndex.load(splittingIndex.toPath());
        Assert.assertTrue(splittingBAMIndex.size() > 0 );
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testNegativeGranularity(){
        final ArgumentsBuilder args = getInputAndOutputArgs(SORTED_BAM, getTempIndexFile())
            .add("--"+ CreateHadoopBamSplittingIndex.SPLITTING_INDEX_GRANULARITY_LONG_NAME).add(-10);
        this.runCommandLine(args);
    }

    @DataProvider(name="unindexable")
    public Object[][] getUnindexableFiles(){
        return new Object[][]{
                {getTestFile("count_reads.cram")},
                {getTestFile("count_reads.sam")}
        };
    }

    @Test(expectedExceptions = UserException.BadInput.class, dataProvider = "unindexable")
    public void testUnindexableFilesFail(final File badFile){
        final ArgumentsBuilder args = getInputAndOutputArgs(badFile, getTempIndexFile());
        this.runCommandLine(args);
    }

    @Test(dataProvider = "indexable")
    public void testUnspecifiedOutputProducesAdjacentIndex(final File bam) throws IOException {
        // copy the bam file to a new temp file
        // we're going to write an index next to it on disk, and we don't want to write into the test resources folder
        final File bamCopy = createTempFile("copy-"+bam, ".bam");
        Files.copy(bam.toPath(), bamCopy.toPath(), StandardCopyOption.REPLACE_EXISTING);
        final File expectedIndex = new File(bamCopy.toPath() + SBIIndex.FILE_EXTENSION);
        Assert.assertFalse(expectedIndex.exists());
        final ArgumentsBuilder args = new ArgumentsBuilder().addInput(bamCopy);
        this.runCommandLine(args);
        expectedIndex.deleteOnExit();
        assertIndexIsNotEmpty(expectedIndex);
    }

    @Test
    public void testBothPathsProduceSameIndex() throws IOException {
        final File splittingIndexOnly = getTempIndexFile();
        final ArgumentsBuilder splittingIndexOnlyArgs = getInputAndOutputArgs(SORTED_BAM, splittingIndexOnly);
        this.runCommandLine(splittingIndexOnlyArgs);
        assertIndexIsNotEmpty(splittingIndexOnly);

        final File splittingIndexWithBai = getTempIndexFile();
        final ArgumentsBuilder splittingAndBaiArgs = getInputAndOutputArgs(SORTED_BAM, splittingIndexWithBai)
                .add("--"+ CreateHadoopBamSplittingIndex.CREATE_BAI_LONG_NAME);
        this.runCommandLine(splittingAndBaiArgs);
        assertIndexIsNotEmpty(splittingIndexWithBai);

        IOUtil.assertFilesEqual(splittingIndexOnly, splittingIndexWithBai);
    }

    public ArgumentsBuilder getInputAndOutputArgs(final File inputFile, final File splittingIndexWithBai) {
        return new ArgumentsBuilder()
                .addInput(inputFile)
                .addOutput(splittingIndexWithBai);
    }

    @Test
    public void testCreateWithBaiCreatesBai(){
        final File splittingIndex = getTempIndexFile();
        final File baiIndex = IOUtils.replaceExtension(splittingIndex, BAMIndex.BAMIndexSuffix);
        Assert.assertFalse(baiIndex.exists());
        final ArgumentsBuilder args = getInputAndOutputArgs(SORTED_BAM, splittingIndex)
                .add("--" + CreateHadoopBamSplittingIndex.CREATE_BAI_LONG_NAME);
        this.runCommandLine(args);
        Assert.assertTrue(baiIndex.exists());
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testCantCreateBaiForUnsortedFile(){
        final ArgumentsBuilder args = getInputAndOutputArgs(UNSORTED_BAM, getTempIndexFile())
                .add("--"+ CreateHadoopBamSplittingIndex.CREATE_BAI_LONG_NAME);
        this.runCommandLine(args);
    }

    private static File getTempIndexFile() {
        return createTempFile("index", "bam" + SBIIndex.FILE_EXTENSION);
    }


}
