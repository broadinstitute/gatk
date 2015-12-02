package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.util.FastqQualityFormat;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Tests for FastqToBam
 */
public final class FastqToSamIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/conversion/fastq2bam");

    public String getTestedClassName() {
        return FastqToSam.class.getSimpleName();
    }

    // fastq files with legal values for each fastq version
    @DataProvider(name = "okVersionFiles")
    public Object[][] okVersionFiles() {
        return new Object[][] {
                {"fastq-sanger/5k-v1-Rhodobacter_LW1.sam.fastq",      FastqQualityFormat.Standard },
                {"fastq-sanger/5k-30BB2AAXX.3.aligned.sam.fastq",     FastqQualityFormat.Standard },
                {"fastq-sanger/sanger_full_range_as_sanger-63.fastq", FastqQualityFormat.Standard }, // all sanger chars

                {"fastq-solexa/s_1_sequence.txt", FastqQualityFormat.Solexa},
                {"fastq-solexa/solexa_full_range_as_solexa.fastq", FastqQualityFormat.Solexa}, // all solexa chars

                {"fastq-illumina/s_1_2_sequence.txt", FastqQualityFormat.Illumina},
                {"fastq-illumina/s_1_2_sequence.txt", FastqQualityFormat.Illumina},
                {"fastq-illumina/s_1_sequence.txt", FastqQualityFormat.Illumina},
        };
    }

    // Illegal values fastq files for each fastq version
    @DataProvider(name = "badVersionFiles")
    public Object[][] badVersionFiles() {
        return new Object[][] {
                {"fastq-sanger/sanger_full_range_as_sanger-63.fastq", FastqQualityFormat.Illumina},
                {"fastq-solexa/s_1_sequence.txt", FastqQualityFormat.Illumina},
        };
    }

    // Illegal fastq format, i.e. doesn't contain correct four lines per record
    @DataProvider(name = "badFormatFiles")
    public Object[][] badFormatFiles() {
        return new Object[][] {
                {"bad-format/bad-qual-header.txt"},
                {"bad-format/bad-seq-header.txt"},
                {"bad-format/extra-line.txt"},
                {"bad-format/too-many-quals.txt"},
                {"bad-format/1lines.txt"},
                {"bad-format/2lines.txt"},
                {"bad-format/3lines.txt"},
        };
    }

    // permissive fastq format, i.e. contain blank lines at various places
    @DataProvider(name = "permissiveFormatFiles")
    public Object[][] permissiveFormatFiles() {
        return new Object[][] {
                {"permissive-format/pair1.txt",          "permissive-format/pair2.txt", FastqQualityFormat.Standard },
                {"permissive-format/s_1_1_sequence.txt",    "permissive-format/s_1_2_sequence.txt", FastqQualityFormat.Illumina},
                {"permissive-format/pair1.txt", null, FastqQualityFormat.Standard},
                {"permissive-format/pair2.txt", null, FastqQualityFormat.Standard},
                {"permissive-format/s_1_1_sequence.txt", null, FastqQualityFormat.Illumina},
                {"permissive-format/s_1_2_sequence.txt", null, FastqQualityFormat.Illumina},
                {"permissive-format/s_1_sequence.txt", null, FastqQualityFormat.Illumina},

        };
    }


    // OK paired fastq files
    @DataProvider(name = "okPairedFiles")
    public Object[][] okPairedFiles() {
        return new Object[][] {
                {"ok-paired/pair1.txt",          "ok-paired/pair2.txt", FastqQualityFormat.Standard },
                {"fastq-illumina/s_1_1_sequence.txt", "fastq-illumina/s_1_2_sequence.txt", FastqQualityFormat.Illumina}
        };
    }

    // Inconsistent paired fastq files
    @DataProvider(name = "badPairedFiles")
    public Object[][] badPairedFiles() {
        return new Object[][] {
                {"ok-paired/pair1.txt",          "bad-paired/pair2-one-more-record.txt" },
                {"bad-paired/pair1-one-more-record.txt", "ok-paired/pair2.txt" },
                {"ok-paired/pair1.txt",          "bad-paired/pair2-badnum.txt" },
                {"bad-paired/pair1-badnum.txt",   "ok-paired/pair2.txt" },
                {"bad-paired/pair1-nonum.txt",    "ok-paired/pair2.txt" },
                {"bad-paired/pair1-onetoken.txt", "ok-paired/pair2.txt" },
        };
    }

    @Test(dataProvider = "permissiveFormatFiles")
    public void testPermissiveOk(final String filename1, final String filename2, final FastqQualityFormat version) throws IOException {
        convertFile(filename1,filename2,version,true);
    }

    @Test(dataProvider = "permissiveFormatFiles",expectedExceptions = SAMException.class)
    public void testPermissiveFail(final String filename1, final String filename2, final FastqQualityFormat version) throws IOException {
        convertFile(filename1,filename2,version,false);
    }

    @Test(dataProvider = "okVersionFiles")
    public void testFastqVersionOk(final String fastqVersionFilename, final FastqQualityFormat version) throws IOException {
        convertFile(fastqVersionFilename, version);
    }

    @Test(dataProvider = "badVersionFiles", expectedExceptions = SAMException.class)
    public void testFastqVersionBad(final String fastqVersionFilename, final FastqQualityFormat version) throws IOException {
        convertFile(fastqVersionFilename, version);
    }

    @Test(dataProvider = "badFormatFiles", expectedExceptions = SAMException.class)
    public void testBadFile(final String filename) throws IOException {
        convertFile(filename, null, FastqQualityFormat.Standard);
    }

    @Test(dataProvider = "badPairedFiles", expectedExceptions = UserException.class)
    public void testPairedBad(final String filename1, final String filename2) throws IOException {
        convertFile(filename1, filename2, FastqQualityFormat.Standard);
    }

    @Test(dataProvider = "okPairedFiles")
    public void testPairedOk(final String filename1, final String filename2, final FastqQualityFormat version) throws IOException {
        convertFile(filename1, filename2, version);
    }

    private File convertFile(final String filename, final FastqQualityFormat version) throws IOException {
        return convertFile(filename, null, version);
    }

    private File convertFile(final String fastqFilename1, final String fastqFilename2, final FastqQualityFormat version) throws IOException{
        return convertFile(fastqFilename1, fastqFilename2, version,false);
    }

    private File convertFile(final String fastqFilename1, final String fastqFilename2, final FastqQualityFormat version,final boolean permissiveFormat) throws IOException {
        final File fastq1 = new File(TEST_DATA_DIR, fastqFilename1);
        final File fastq2 = (fastqFilename2 != null) ? new File(TEST_DATA_DIR, fastqFilename2) : null;
        final File samFile = newTempSamFile(fastq1.getName());

        final List<String> args = new ArrayList<>();

        args.add("--FASTQ");
        args.add(fastq1.getAbsolutePath());
        if (fastqFilename2 != null) {
            args.add("--FASTQ2");
            args.add(fastq2.getAbsolutePath());
        }
        args.add("--output");
        args.add(samFile.getAbsolutePath());
        args.add("--QUALITY_FORMAT");
        args.add(version.toString());
        args.add("--READ_GROUP_NAME");
        args.add("rg");
        args.add("--SAMPLE_NAME");
        args.add("s1");
        if (permissiveFormat) {
            args.add("--ALLOW_AND_IGNORE_EMPTY_LINES");
            args.add("true");
        }

        runCommandLine(args);
        return samFile;
    }

    private static File newTempSamFile(final String filename) throws IOException {
        final File file = BaseTest.createTempFile(filename, ".sam");
        return file;
    }

    private static File newTempFile(final String filename) throws IOException {
        final File file = BaseTest.createTempFile(filename, ".tmp");
        return file;
    }

    //  Test for legal syntax for pair read names for FastqToSam.getBaseName()
    //  We create a dummy file to test the getBaseName() method since it expects
    //  an existing file.

    // TODO - Should switch over to using invocation via new PicardCommandLine() - BUT the tests using this are accessing class members directly.
    private static final FastqToSam fastqToSam = new FastqToSam();
    private static FastqReader freader1;
    private static FastqReader freader2;

    @BeforeClass
    public static void beforeClass() throws IOException {
        final File dummyFile = newTempFile("dummy");
        freader1 = new FastqReader(dummyFile);
        freader2 = new FastqReader(dummyFile);
    }

    @DataProvider(name = "okPairNames")
    public Object[][] okPairNames() {
        return new Object[][] {
                {"aa/1", "aa/2" },
                {"aa", "aa" },
                {"aa/bb", "aa/bb" },
                {"aa/bb/", "aa/bb/" },
                {"aa/bb/1", "aa/bb/2" },
                {"aa/bb/cc/dd/ee/ff/1", "aa/bb/cc/dd/ee/ff/2" },
                {"////1", "////2" },
                {"/", "/" },
                {"////", "////" },
                {"/aa", "/aa" },
                {"aa/", "aa/" },
                {"ab/c", "ab/c"}
        };
    }

    @DataProvider(name = "badPairNames")
    public Object[][] badPairNames() {
        return new Object[][] {
                {"", "" },
                {"aa/1", "bb/2" },
                {"aa"  , "bb" },
                {"aa/1", "aa" },
                {"aa",   "aa/2" },
                {"aa/1", "aa/1" },
                {"aa/2", "aa/2" },
        };
    }

    @Test(dataProvider = "okPairNames")
    public void readPairNameOk(final String name1, final String name2) throws IOException {
        fastqToSam.getBaseName(name1, name2, freader1, freader2);
    }

    @Test(dataProvider = "badPairNames", expectedExceptions = {UserException.class, GATKException.class})
    public void readPairNameBad(final String name1, final String name2) throws IOException {
        fastqToSam.getBaseName(name1, name2, freader1, freader2);
    }
}
