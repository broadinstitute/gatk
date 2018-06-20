package org.broadinstitute.hellbender.engine;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.PrintReads;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class CRAMSupportIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/engine");

    @Override
    public String getTestedClassName() {
        return PrintReads.class.getSimpleName();
    }

    @DataProvider(name = "ReadEntireCramTestData")
    public Object[][] readEntireCramTestData() {
        final File ref = new File(hg19MiniReference);
        final List<Object[]> testCases = new ArrayList<>();
        for ( final String outputExtension : Arrays.asList(".sam", ".bam", ".cram") ) {
            // cram, reference, output file extension, expected read names
            testCases.add(new Object[]{new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k")});
            testCases.add(new Object[]{new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k")});
        }
        return testCases.toArray(new Object[][]{});
    }

    @DataProvider(name = "cramFromPathIndex")
    public Object[][] cramFromPathIndex() {
        final File ref = new File(hg19MiniReference);
        return new Object[][] {
            // cram, index file, reference, intervals, expected read names
            { new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), new File(TEST_DATA_DIR, "cram_with_bai_index.cram.bai"),
                    ref, Arrays.asList("1:1-200", "1:1000-2000"), Arrays.asList("a", "d", "e") },
            { new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), new File(TEST_DATA_DIR, "cram_with_crai_index.cram.crai"),
                    ref, Arrays.asList("1:1000-2000", "3:300-400"), Arrays.asList("d", "e", "i", "j") }
        };
    }

    @Test(dataProvider = "cramFromPathIndex")
    public void testQueryCRAMWithIndexFromPath(
            final File cramFile,
            final File indexFile,
            final File referenceFile,
            final List<String> intervalArgs,
            final List<String> expectedReadNames ) throws IOException {
        final List<String> args = new ArrayList<>();
        try (final FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path jimfsRootPath = jimfs.getRootDirectories().iterator().next();

            // copy the cram, the index, the reference, the reference dictionary and the reference index to jimfs
            final Path cramPath = Files.copy(cramFile.toPath(), jimfsRootPath.resolve(cramFile.getName()));
            Files.copy(indexFile.toPath(), jimfsRootPath.resolve(indexFile.getName()));

            final Path refPath = Files.copy(referenceFile.toPath(), jimfsRootPath.resolve(referenceFile.getName()));

            final File referenceDict = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(referenceFile);
            final Path referenceDictPath = referenceDict.toPath();
            Files.copy(referenceDictPath, jimfsRootPath.resolve(referenceDict.getName()));

            final Path referenceIndexPath = ReferenceSequenceFileFactory.getFastaIndexFileName(referenceFile.toPath());
            Files.copy(referenceIndexPath, jimfsRootPath.resolve(referenceFile.getName() + ".fai"));

            final File outputFile = createTempFile("testReadCramWithIndexFromPath", ".bam");
            args.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
            args.add(cramPath.toAbsolutePath().toUri().toString());
            args.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
            args.add(outputFile.getAbsolutePath());
            args.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
            args.add(refPath.toAbsolutePath().toUri().toString());
            intervalArgs.stream().forEach(intervalArg -> { args.add("-L"); args.add(intervalArg); });

            runCommandLine(args);

            SamAssertionUtils.assertCRAMContentsIfCRAM(outputFile);
            checkReadNames(outputFile.toPath(), refPath, expectedReadNames);
        }

    }

    @Test(dataProvider = "ReadEntireCramTestData")
    public void testReadEntireCram(final File cramFile, final File reference, final String outputExtension, final List<String> expectedReadNames ) throws IOException {
        final File outputFile = createTempFile("testReadEntireCram", outputExtension);
        final String[] args = new String[] {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, cramFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, reference.getAbsolutePath()
        };
        runCommandLine(args);

        SamAssertionUtils.assertCRAMContentsIfCRAM(outputFile);
        checkReadNames(outputFile.toPath(), reference.toPath(), expectedReadNames);
    }

    @Test
    public void testSamtoolsGeneratedCRAMSliceMD5Calculation() throws IOException {
        // Note: The input CRAM used for this test was generated using samtools, with the "ambiguityCodes.fasta" file
        // as the reference. Since the reference contains ambiguity codes (essential, since part of what we're trying
        // to validate is that htsjdk calculates the same MD5 value for a slice who's reference spans ambiguity codes
        // as samtools does), and since GATK wants an accompanying sequence dictionary, the .dict was generated
        // with GATK because samtools has a bug in dictionary generation when the reference has ambiguity codes.
        // See https://github.com/samtools/samtools/issues/704, and gatk tracking issue:
        // https://github.com/broadinstitute/gatk/issues/3306
        final File samtoolsGeneratedCRAM = new File(TEST_DATA_DIR, "samtoolsSliceMD5WithAmbiguityCodesTest.cram");
        final File referenceWithAmbiguityCodes = new File(TEST_DATA_DIR, "ambiguityCodes.fasta");
        final File outputFile = createTempFile("testReadSamtoolsGeneratedCRAM", ".cram");
        final String[] args = new String[] {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, samtoolsGeneratedCRAM.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, referenceWithAmbiguityCodes.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.DISABLE_SEQUENCE_DICT_VALIDATION_NAME
        };
        runCommandLine(args);  // no assert, just make sure we don't throw
    }


    @DataProvider(name = "ReadCramWithIntervalsIndexTestData")
    public Object[][] readCramWithIntervalsBAIIndexTestData() {
        final File ref = new File(hg19MiniReference);
        final List<Object[]> testCases = new ArrayList<>();
        for ( final String outputExtension : Arrays.asList(".sam", ".bam", ".cram") ) {
            // cram, reference, output file extension, intervals, expected read names
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000"), Arrays.asList("d", "e") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-1099"), Arrays.asList("d") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("1:1-200", "1:1000-2000"), Arrays.asList("a", "d", "e") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000", "3:300-400"), Arrays.asList("d", "e", "i", "j") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_bai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000", "3:300-399"), Arrays.asList("d", "e", "i") });

            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000"), Arrays.asList("d", "e") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-1099"), Arrays.asList("d") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("1:1-200", "1:1000-2000"), Arrays.asList("a", "d", "e") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000", "3:300-400"), Arrays.asList("d", "e", "i", "j") });
            testCases.add(new Object[]{ new File(TEST_DATA_DIR, "cram_with_crai_index.cram"), ref, outputExtension, Arrays.asList("1:1000-2000", "3:300-399"), Arrays.asList("d", "e", "i") });
        }
        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ReadCramWithIntervalsIndexTestData")
    public void testReadCramWithIntervalsWithBAIIndex( final File cramFile, final File reference, final String outputExtension,
                                                       final List<String> intervalArgs, final List<String> expectedReadNames ) throws IOException {
        final File outputFile = createTempFile("testReadCramWithIntervalsWithBAIIndex", outputExtension);
        final List<String> args = new ArrayList<>();
        args.addAll(Arrays.asList(
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, cramFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, reference.getAbsolutePath()
        ));
        intervalArgs.stream().forEach(intervalArg -> { args.add("-L"); args.add(intervalArg); });

        runCommandLine(args);

        SamAssertionUtils.assertCRAMContentsIfCRAM(outputFile);
        checkReadNames(outputFile.toPath(), reference.toPath(), expectedReadNames);
    }

    private void checkReadNames( final Path outputPath, final Path reference, final List<String> expectedReadNames ) throws IOException {
        List<String> actualReadNames = new ArrayList<>();
        try ( final SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).referenceSequence(reference).open(outputPath) ) {
            for ( SAMRecord read : reader ) {
                actualReadNames.add(read.getReadName());
            }
        }
        Assert.assertEquals(actualReadNames, expectedReadNames, "Read names in output do not match expected read names");
    }

    @Test(dataProvider="testingDataNoRef", expectedExceptions = UserException.MissingReference.class)
    public void testNoRef(String fileIn, String extOut) throws Exception {
        final File outFile = GATKBaseTest.createTempFile(fileIn + ".", extOut);
        File readInput = new File(TEST_DATA_DIR, fileIn);
        final String[] args = new String[]{
                "--input" , readInput.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);
    }

    @DataProvider(name="testingDataNoRef")
    public Object[][] testingDataNoRef() {
        return new String[][]{
                {"cramtest.cram", ".sam"}
        };
    }

    @Test(dataProvider="testingDataNoRefMultipleInputs", expectedExceptions = UserException.MissingReference.class)
    public void testNoRefMulti(String fileIn1, String fileIn2, String extOut) throws Exception {
        final File outFile = GATKBaseTest.createTempFile(fileIn1 + ".", extOut);
        File readInput1 = new File(TEST_DATA_DIR, fileIn1);
        File readInput2 = new File(TEST_DATA_DIR, fileIn2);
        final String[] args = new String[]{
                "--input" , readInput1.getAbsolutePath(),
                "--input" , readInput2.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);
    }

    @DataProvider(name="testingDataNoRefMultipleInputs")
    public Object[][] testingDataNoRefMulti() {
        return new String[][]{
                {"cramtest.sam", "cramtest.cram", ".bam"}
        };
    }

    // This test case shows that when a CRAM input is provided with a reference that does not have all of the contigs
    // from the CRAM in its sequence dictionary, we throw a UserException.
    @Test(dataProvider="testingDataWrongRef", expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testWrongRef(String fileIn, String extOut, String referenceFile) throws Exception {
        final File outFile = GATKBaseTest.createTempFile(fileIn + ".", extOut);
        File readInput = new File(TEST_DATA_DIR, fileIn);
        File reference = new File(TEST_DATA_DIR, referenceFile);
        final String[] args = new String[]{
                "--input" , readInput.getAbsolutePath(),
                "--output", outFile.getAbsolutePath(),
                "-R", reference.getAbsolutePath()
        };
        runCommandLine(args);
    }

    @DataProvider(name="testingDataWrongRef")
    public Object[][] testingDataWrongRef() {
        return new String[][]{
                {"cramtest.cram", ".sam", "cramtestWrongRef.fasta"},
        };
    }


}