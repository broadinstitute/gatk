package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import org.apache.commons.io.FileUtils;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadLengthReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadNameReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public final class PrintReadsIntegrationTest extends CommandLineProgramTest{

    private static final File TEST_DATA_DIR = getTestDataDir();

    @Override
    public String getTestedClassName() {
        return PrintReads.class.getSimpleName();
    }

    public void doFileToFile(String fileIn, String extOut, String reference, boolean testMD5) throws Exception {
        String samFile = fileIn;
        final File outFile = BaseTest.createTempFile(samFile + ".", extOut);
        final File ORIG_BAM = new File(TEST_DATA_DIR, samFile);
        final File refFile;

        final ArrayList<String> args = new ArrayList<>();
        args.add("--input"); args.add(ORIG_BAM.getAbsolutePath());
        args.add("--output"); args.add(outFile.getAbsolutePath());
        if (reference != null) {
            refFile = new File(TEST_DATA_DIR, reference);
            args.add("-R"); args.add(refFile.getAbsolutePath());
        }
        else {
            refFile = null;
        }
        if (testMD5) {
            args.add("--createOutputBamMD5");
            args.add("true");
        }
        runCommandLine(args);

        SamAssertionUtils.assertSamsEqual(outFile, ORIG_BAM, refFile);

        if (testMD5) {
            checkMD5asExpected(outFile);
        }
    }

    private void checkMD5asExpected(final File outFile) throws IOException {
        final File md5File = new File(outFile.getAbsolutePath() + ".md5");
        if (md5File.exists()) {
            md5File.deleteOnExit();
        }
        Assert.assertTrue(md5File.exists(), md5File + " does not exist");
        final String expectedMD5 = Utils.calculateFileMD5(outFile);
        final String actualMD5 = FileUtils.readFileToString(md5File);
        Assert.assertEquals(actualMD5, expectedMD5);
    }

    @Test(dataProvider="testingData")
    public void testFileToFile(String fileIn, String extOut, String reference) throws Exception {
        doFileToFile(fileIn, extOut, reference, false);
    }

    @Test(dataProvider="testingData")
    public void testFileToFileWithMD5(String fileIn, String extOut, String reference) throws Exception {
        doFileToFile(fileIn, extOut, reference, true);
    }

    @DataProvider(name="testingData")
    public Object[][] testingData() {
        return new String[][]{
                {"print_reads.sam", ".sam", null},
                {"print_reads.sam", ".bam", null},
                {"print_reads.sam", ".cram", "print_reads.fasta"},
                {"print_reads.bam", ".sam", null},
                {"print_reads.bam", ".bam", null},
                {"print_reads.bam", ".cram", "print_reads.fasta"},
                {"print_reads.cram", ".sam", "print_reads.fasta"},
                {"print_reads.cram", ".bam", "print_reads.fasta"},
                {"print_reads.cram", ".cram", "print_reads.fasta"},

                {"print_reads.sorted.sam", ".sam", null},
                {"print_reads.sorted.sam", ".bam", null},
                {"print_reads.sorted.sam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.bam", ".sam", null},
                {"print_reads.sorted.bam", ".bam", null},
                {"print_reads.sorted.bam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".sam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".bam", "print_reads.fasta"},
                {"print_reads.sorted.cram", ".cram", "print_reads.fasta"},

                {"print_reads.sorted.queryname.sam", ".sam", null},
                {"print_reads.sorted.queryname.sam", ".bam", null},
                {"print_reads.sorted.queryname.sam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.bam", ".sam", null},
                {"print_reads.sorted.queryname.bam", ".bam", null},
                {"print_reads.sorted.queryname.bam", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".sam", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".bam", "print_reads.fasta"},
                {"print_reads.sorted.queryname.cram", ".cram", "print_reads.fasta"},

                //test queryname-sorted crams with multiref containers in GATK:
                //print_reads.sorted.queryname_htsjdk_2.1.0.cram was generated from print_reads.sam
                //using gatk4 PrintReads/htsjdk.2.1.0, which includes changes to support
                //multireference containers
                {"print_reads.sorted.queryname.htsjdk-2.1.0.cram", ".cram", "print_reads.fasta"},
                {"print_reads.sorted.queryname.htsjdk-2.1.0.cram", ".sam", "print_reads.fasta"}
        };
    }

    @Test
    public void testReadThatConsumesNoReferenceBases() throws IOException {
        final File zeroRefBasesReadBam = new File(TEST_DATA_DIR, "read_consumes_zero_ref_bases.bam");
        final File outFile = BaseTest.createTempFile("testReadThatConsumesNoReferenceBases", ".bam");
        final String[] args = new String[] {
                "--input" , zeroRefBasesReadBam.getAbsolutePath(),
                "--output", outFile.getAbsolutePath()
        };
        // Make sure no exception is thrown given an input containing a read that consumes no reference bases
        runCommandLine(args);

        //Make sure we print the read, ie not lose it.
        SamAssertionUtils.assertSamsEqual(outFile, zeroRefBasesReadBam);
    }

    @Test
    public void testNoConflictPG() throws IOException {
        final File inFile = new File(TEST_DATA_DIR, "print_reads_withPG.sam");
        final File outFile = BaseTest.createTempFile("testNoConflictRG", ".sam");
        final String[] args = new String[] {
                "--input" , inFile.getAbsolutePath(),
                "--addOutputSAMProgramRecord",
                "--output", outFile.getAbsolutePath()
        };
        runCommandLine(args);

        //Make sure contents are the same
        SamAssertionUtils.assertSamsEqual(outFile, inFile);

        //input has GATK PrintReads not NOT GATK PrintReads.1 in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK PrintReads"));
        Assert.assertNull(SamReaderFactory.makeDefault().open(inFile).getFileHeader().getProgramRecord("GATK PrintReads.1"));

        //output has both GATK PrintReads and GATK PrintReads.1 in headers
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK PrintReads"));
        Assert.assertNotNull(SamReaderFactory.makeDefault().open(outFile).getFileHeader().getProgramRecord("GATK PrintReads.1"));
    }

    @DataProvider(name = "UnmappedReadInclusionTestData")
    public Object[][] unmappedReadInclusionTestData() {
        // This bam has mapped reads from various contigs, plus a few unmapped reads with no mapped mate
        final File unmappedBam = new File(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_with_unmapped.bam");

        // This is a snippet of the CEUTrio.HiSeq.WGS.b37.NA12878 bam from large, with mapped reads
        // from chromosome 20 (with one mapped read having an unmapped mate), plus several unmapped
        // reads with no mapped mate.
        final File ceuSnippet = new File(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.bam");
        final File ceuSnippetCram = new File(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.cram");

        return new Object[][] {
                { unmappedBam, null, Arrays.asList("unmapped"), Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                // The same interval as above in an intervals file
                { unmappedBam, null, Arrays.asList(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_unmapped.intervals"), Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                { unmappedBam, null, Arrays.asList("1:200-300", "unmapped"), Arrays.asList("a", "b", "c", "u1", "u2", "u3", "u4", "u5") },
                { unmappedBam, null, Arrays.asList("1:200-300", "4:700-701", "unmapped"), Arrays.asList("a", "b", "c", "k", "u1", "u2", "u3", "u4", "u5") },
                // The same intervals as above in an intervals file
                { unmappedBam, null, Arrays.asList(TestResources.publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_unmapped2.intervals"), Arrays.asList("a", "b", "c", "k", "u1", "u2", "u3", "u4", "u5") },
                { ceuSnippet, null, Arrays.asList("unmapped"), Arrays.asList("g", "h", "h", "i", "i") },
                { ceuSnippet, null, Arrays.asList("20:10000009-10000011", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "g", "h", "h", "i", "i") },
                { ceuSnippet, null, Arrays.asList("20:10000009-10000013", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "f", "f", "g", "h", "h", "i", "i") },
                { ceuSnippetCram, TestResources.b37_reference_20_21, Arrays.asList("unmapped"), Arrays.asList("g", "h", "h", "i", "i") },
                { ceuSnippetCram, TestResources.b37_reference_20_21, Arrays.asList("20:10000009-10000011", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "g", "h", "h", "i", "i") },
                { ceuSnippetCram, TestResources.b37_reference_20_21, Arrays.asList("20:10000009-10000013", "unmapped"), Arrays.asList("a", "b", "c", "d", "e", "f", "f", "g", "h", "h", "i", "i") }
        };
    }

    @Test(dataProvider = "UnmappedReadInclusionTestData")
    public void testUnmappedReadInclusion( final File input, final String reference, final List<String> intervalStrings, final List<String> expectedReadNames ) {
        final File outFile = createTempFile("testUnmappedReadInclusion", ".bam");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("-I"); args.add(input.getAbsolutePath());
        args.add("-O"); args.add(outFile.getAbsolutePath());
        for ( final String intervalString : intervalStrings ) {
            args.add("-L"); args.add(intervalString);
        }
        if ( reference != null ) {
            args.add("-R"); args.add(reference);
        }

        runCommandLine(args);

        try ( final ReadsDataSource outputReadsSource = new ReadsDataSource(outFile.toPath()) ) {
            final List<GATKRead> actualReads = new ArrayList<>();
            for ( final GATKRead read : outputReadsSource ) {
                actualReads.add(read);
            }

            Assert.assertEquals(actualReads.size(), expectedReadNames.size(), "Wrong number of reads output");

            for ( int readNumber = 0; readNumber < actualReads.size(); ++readNumber ) {
                Assert.assertEquals(actualReads.get(readNumber).getName(), expectedReadNames.get(readNumber), "Unexpected read name");
            }
        }
    }

    @DataProvider(name="readFilterTestData")
    public Object[][] testReadFilterData() {
        return new Object[][]{
                {"print_reads_one_malformed_read.sam", null, ".sam", Collections.emptyList(), 7},
                {"print_reads_one_malformed_read.sam", null, ".sam", Arrays.asList("--disableToolDefaultReadFilters"), 8},
                {"print_reads_one_malformed_read.sam", null, ".sam",
                        Arrays.asList("--disableReadFilter", "WellformedReadFilter"), 8},
                {"print_reads.sorted.sam", null, ".sam", Arrays.asList("--disableToolDefaultReadFilters"), 8},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--readFilter", ReadNameReadFilter.class.getSimpleName(),
                                "--readName", "both_reads_align_clip_adapter"),
                        2},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "100",
                                "--maxReadLength", "200"),
                        8},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "1",
                                "--maxReadLength", "10"),
                        0},
                {"print_reads.sorted.sam", null, ".sam",
                        Arrays.asList(
                                "--readFilter", ReadNameReadFilter.class.getSimpleName(),
                                "--readName", "both_reads_align_clip_adapter",
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "100",
                                "--maxReadLength", "101"),
                        2},
                {"print_reads.sorted.bam", null, ".sam", Arrays.asList("--disableToolDefaultReadFilters"), 8},
                {"print_reads.sorted.bam", null, ".sam",
                        Arrays.asList(
                                "--readFilter", ReadNameReadFilter.class.getSimpleName(),
                                "--readName", "both_reads_align_clip_adapter"),
                        2},
                {"print_reads.sorted.bam", null, ".sam",
                        Arrays.asList(
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "100",
                                "--maxReadLength", "101"),
                        8},
                {"print_reads.sorted.bam", null, ".sam",
                        Arrays.asList(
                                "--readFilter", ReadNameReadFilter.class.getSimpleName(),
                                "--readName", "both_reads_align_clip_adapter",
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "100",
                                "--maxReadLength", "101"),
                        2},
                {"print_reads.sorted.cram", "print_reads.fasta", ".sam",
                        Arrays.asList(
                                "--readFilter", ReadNameReadFilter.class.getSimpleName(),
                                "--readName", "both_reads_align_clip_adapter",
                                "--RF", ReadLengthReadFilter.class.getSimpleName(),
                                "--minReadLength", "100",
                                "--maxReadLength", "101"),
                        2},
        };
    }

    @Test(dataProvider = "readFilterTestData")
    public void testReadFilters(
            final String input,
            final String reference,
            final String extOut,
            final List<String> inputArgs,
            final int expectedCount) throws IOException
    {
        final File outFile = createTempFile("testReadFilter", extOut);

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("-I"); args.add(new File(TEST_DATA_DIR, input).getAbsolutePath());
        args.add("-O"); args.add(outFile.getAbsolutePath());
        if ( reference != null ) {
            args.add("-R"); args.add(new File(TEST_DATA_DIR, reference).getAbsolutePath());
        }
        for (final String filter : inputArgs) {
            args.add(filter);
        }

        runCommandLine(args);


        SamReaderFactory factory = SamReaderFactory.makeDefault();
        if (reference != null) {
            factory = factory.referenceSequence(new File(TEST_DATA_DIR, reference));
        }
        int count = 0;
        try (final SamReader reader = factory.open(outFile)) {
            Iterator<SAMRecord> it = reader.iterator();
            while (it.hasNext()) {
                SAMRecord rec = it.next();
                count++;
            }
        }
        Assert.assertEquals(count, expectedCount);
    }

    @Test(expectedExceptions = UserException.MissingReference.class)
    public void testNonExistentReference() throws Exception {
        final File inCram = new File(TEST_DATA_DIR, "print_reads.sorted.cram");
        final File outCram = BaseTest.createTempFile("print_reads_bad_reference", ".cram");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(inCram.getCanonicalPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
        args.add(outCram.getCanonicalPath());
        args.add("-R");
        args.add(BaseTest.getSafeNonExistentFile("Nonexistent.fasta").getCanonicalPath());

        runCommandLine(args.getArgsArray());
    }
}