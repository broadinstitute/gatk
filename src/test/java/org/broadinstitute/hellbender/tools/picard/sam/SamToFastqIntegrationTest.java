package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * Tests for SamToFastq
 */
public final class SamToFastqIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/conversion/bam2fastq/paired");
    private static final String CLIPPING_TEST_DATA = "ok/clipping_test.sam";

    public String getTestedClassName() {
        return SamToFastq.class.getSimpleName();
    }

    @DataProvider(name = "okFiles")
    public Object[][] okFiles() {
        return new Object[][] {
                {"ok/sorted-pair.sam"}, // 5 sorted pairs (10 records) - mate1, mate2
                {"ok/sorted-pair-no-rg.sam"}, // 5 sorted pairs (10 records) - mate1, mate2, no read group
                {"ok/last-pair-mates-flipped.sam" }, // 4 pairs, :05 mate2, mate1
                {"ok/first-mate-bof-last-mate-eof.sam"}, // :01 mate1, 4 pairs, :01 mate2
        };
    }


    @DataProvider(name = "badFiles")
    public Object[][] badFiles() {
        return new Object[][] {
                {"bad/unpaired-mate.sam"} // mate1 without its mate2
        };
    }

    private void convertFile(final String [] args) {
        runCommandLine(args);
    }

    @Test(dataProvider = "clippingTests")
    public void testClipping(final String clippingAction, final String bases1_1, final String quals1_1, final String bases1_2, final String quals1_2,
                             final String bases2_1, final String quals2_1, final String bases2_2, final String quals2_2, final String testName) throws IOException {
        final File samFile = new File(TEST_DATA_DIR, CLIPPING_TEST_DATA) ;
        final File f1 = BaseTest.createTempFile("clippingtest1", "fastq");
        final File f2 = BaseTest.createTempFile("clippingtest2", "fastq");

        if (clippingAction != null) {
            convertFile(new String[]{
                    "-input",             samFile.getAbsolutePath(),
                    "--FASTQ",            f1.getAbsolutePath(),
                    "--SECOND_END_FASTQ", f2.getAbsolutePath(),
                    "--CLIPPING_ACTION",  clippingAction,
                    "--CLIPPING_ATTRIBUTE",  "XT"
            });
        } else {
            convertFile(new String[]{
                    "--input",             samFile.getAbsolutePath(),
                    "--FASTQ",             f1.getAbsolutePath(),
                    "--SECOND_END_FASTQ",  f2.getAbsolutePath(),
            });
        }

        Iterator<FastqRecord> it = new FastqReader(f1).iterator();
        FastqRecord first = it.next();
        Assert.assertEquals(first.getReadString(), bases1_1, testName);
        Assert.assertEquals(first.getBaseQualityString(), quals1_1, testName);
        FastqRecord second = it.next();
        Assert.assertEquals(second.getReadString(), bases1_2, testName);
        Assert.assertEquals(second.getBaseQualityString(), quals1_2, testName);
        it = new FastqReader(f2).iterator();
        first = it.next();
        Assert.assertEquals(first.getReadString(), bases2_1, testName);
        Assert.assertEquals(first.getBaseQualityString(), quals2_1, testName);
        second = it.next();
        Assert.assertEquals(second.getReadString(), bases2_2, testName);
        Assert.assertEquals(second.getBaseQualityString(), quals2_2, testName);
    }

    @DataProvider(name = "clippingTests")
    public Object[][] clippingTests() {
        return new Object[][] {
                {null, "AAAAAAAAAA", "1111111111", "AAAAAAAAAA", "1111111111", "CCCCCCCCCC", "2222222222", "GGGGGGGGGG", "2222222222", "No clipping test"},
                {"X",  "AAAAAAA",    "1111111",    "AAAAAA",     "111111",     "CCCCCCCC",   "22222222",   "GGGGGG",     "222222",     "Cut clipped bases test"},
                {"N",  "AAAAAAANNN", "1111111111", "AAAAAANNNN", "1111111111", "CCCCCCCCNN", "2222222222", "GGGGGGNNNN", "2222222222", "Mask clipped bases test"},
                {"2",  "AAAAAAAAAA", "1111111###", "AAAAAAAAAA", "111111####", "CCCCCCCCCC", "22222222##", "GGGGGGGGGG", "222222####", "Change clipped qualities test"}
        };
    }

    @Test(dataProvider = "okFiles")
    public void testOkFile(final String samFilename) throws IOException {
        final File samFile = new File(TEST_DATA_DIR,samFilename);
        final File pair1File = newTempFastqFile("pair1");
        final File pair2File = newTempFastqFile("pair2");

        convertFile(new String[]{
                "--input", samFile.getAbsolutePath(),
                "--FASTQ", pair1File.getAbsolutePath(),
                "--SECOND_END_FASTQ", pair2File.getAbsolutePath()
        });

        // Check that paired fastq files are same size
        final Set<String> outputHeaderSet1 = createFastqReadHeaderSet(pair1File);
        final Set<String> outputHeaderSet2 = createFastqReadHeaderSet(pair2File);
        Assert.assertEquals(outputHeaderSet1.size(), outputHeaderSet2.size());

        // Create map of mate pairs from SAM records
        final Map<String,MatePair> map = createSamMatePairsMap(samFile) ;
        Assert.assertEquals(map.size(), outputHeaderSet2.size());

        // Ensure that each mate of each pair in SAM file is in the correct fastq pair file
        for (final Map.Entry<String,MatePair> entry : map.entrySet() ) {
            final MatePair mpair = entry.getValue();
            Assert.assertNotNull(mpair.mate1); // ensure we have two mates
            Assert.assertNotNull(mpair.mate2);
            Assert.assertEquals(mpair.mate1.getReadName(),mpair.mate2.getReadName());
            final String readName = mpair.mate1.getReadName() ;
            Assert.assertTrue(outputHeaderSet1.contains(readName+"/1")); // ensure mate is in correct file
            Assert.assertTrue(outputHeaderSet2.contains(readName+"/2"));
        }
    }


    @Test
    public void testHelp() {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        PrintStream ps= new PrintStream(baos);
        System.setErr(ps);
        convertFile(new String[]{
                "--help",
        });
        final String s = baos.toString();
        String[] lines = s.split("\n");
        for (int i = 0; i < lines.length; i++) {
            Assert.assertFalse(lines[i].contains("OUTPUT_PER_RG (OPRG) OUTPUT_PER_RG (OPRG)"), lines[i]); //Same option twice!
        }
    }


    @Test(dataProvider =  "okFiles")
    public void testOkInterleavedFile(final String samFilename) throws IOException {
        final File samFile = new File(TEST_DATA_DIR,samFilename);
        final File pairFile = newTempFastqFile("pair");

        convertFile(new String[]{
                "--input", samFile.getAbsolutePath(),
                "--FASTQ", pairFile.getAbsolutePath(),
                "--INTERLEAVE",
        });

        final Set<String> outputHeaderSet = createFastqReadHeaderSet(pairFile);
        // Create map of mate pairs from SAM records
        final Map<String,MatePair> map = createSamMatePairsMap(samFile) ;
        Assert.assertEquals(map.size() * 2, outputHeaderSet.size());

        // Ensure that each mate of each pair in SAM file is in the correct fastq pair file
        for (final Map.Entry<String,MatePair> entry : map.entrySet() ) {
            final MatePair mpair = entry.getValue();
            Assert.assertNotNull(mpair.mate1); // ensure we have two mates
            Assert.assertNotNull(mpair.mate2);
            Assert.assertEquals(mpair.mate1.getReadName(),mpair.mate2.getReadName());
            final String readName = mpair.mate1.getReadName() ;
            Assert.assertTrue(outputHeaderSet.contains(readName+"/1")); // ensure mate is in correct file
            Assert.assertTrue(outputHeaderSet.contains(readName+"/2"));
        }
    }


    @Test (dataProvider = "badFiles", expectedExceptions= SAMFormatException.class)
    public void testBadFile(final String samFilename) throws IOException {
        final File samFile = new File(TEST_DATA_DIR,samFilename);
        final File pair1 = BaseTest.createTempFile("tt-pair1.", ".fastq");
        final File pair2 = BaseTest.createTempFile("tt-pair2.", ".fastq");
        convertFile(new String[]{
                "--input", samFile.getAbsolutePath(),
                "--FASTQ", pair1.getAbsolutePath(),
                "--SECOND_END_FASTQ", pair2.getAbsolutePath()
        });
    }

    @DataProvider(name = "okGroupedFiles")
    public Object[][] okGroupedFiles() {
        return new Object[][] {
                {"ok/grouped-last-pair-mates-flipped.sam", null,  null,  new String[]{"rg1","rg2"}},
        };
    }


    @DataProvider(name = "badGroupedFiles")
    public Object[][] badGroupedFiles() {
        return new Object[][] {
                {"bad/grouped-unpaired-mate.sam", null,  null,  new String[]{"rg1.fastq","rg2.fastq"}}
        };
    }

    @Test(dataProvider = "okGroupedFiles")
    public void testOkGroupedFiles(final String samFilename, final String fastq, final String secondEndFastq,
                                   final String [] groupFiles) throws IOException {
        final File samFile = new File(TEST_DATA_DIR,samFilename);
        final Map<String, Set<String>> outputSets = new HashMap<>(groupFiles.length);

        final String tmpDir = IOUtil.getDefaultTmpDir().getAbsolutePath() + "/";
        final String [] args = new String[]{
                "--input", samFile.getAbsolutePath(),
                "--OUTPUT_PER_RG",
                "--OUTPUT_DIR", tmpDir,
        };
        runCommandLine(args);

        File f1;
        File f2;
        String fname1;
        String fname2;
        String keyName1;
        String keyName2;
        Set<String> outputHeaderSet1;
        Set<String> outputHeaderSet2;
        for(final String groupPUName : groupFiles)
        {
            keyName1 = groupPUName + "_1";
            keyName2 = groupPUName + "_2";
            fname1 = tmpDir + "/" + keyName1 + ".fastq";
            fname2 = tmpDir + "/" + keyName2 + ".fastq";
            f1 = new File(fname1);
            f2 = new File(fname2);
            f1.deleteOnExit();
            f2.deleteOnExit();
            IOUtil.assertFileIsReadable(f1);
            IOUtil.assertFileIsReadable(f2);

            // Check that paired fastq files are same size and store them for later comparison
            outputHeaderSet1 = createFastqReadHeaderSet(f1);
            outputHeaderSet2 = createFastqReadHeaderSet(f2);
            outputSets.put(keyName1 , outputHeaderSet1);
            outputSets.put(keyName2, outputHeaderSet2);
            Assert.assertEquals(outputHeaderSet1.size(), outputHeaderSet2.size());
        }

        // Create map of read groups and mate pairs from SAM records
        final Map<String, Map<String,MatePair>> map = createPUPairsMap(samFile);

        for(final Map.Entry<String, Map<String, MatePair>> groupEntry : map.entrySet()) {
            // Ensure that for each group, each mate of each pair in the SAM file is in the correct fastq pair file
            for (final Map.Entry<String,MatePair> entry : groupEntry.getValue().entrySet() ) {
                final MatePair mpair = entry.getValue();
                outputHeaderSet1 = outputSets.get(groupEntry.getKey() + "_1");
                outputHeaderSet2 = outputSets.get(groupEntry.getKey() + "_2");

                Assert.assertNotNull(mpair.mate1); // ensure we have two mates
                Assert.assertNotNull(mpair.mate2);
                Assert.assertEquals(mpair.mate1.getReadName(),mpair.mate2.getReadName());
                final String readName = mpair.mate1.getReadName() ;
                Assert.assertTrue(outputHeaderSet1.contains(readName+"/1")); // ensure mate is in correct file
                Assert.assertTrue(outputHeaderSet2.contains(readName+"/2"));
            }
        }
    }


    @Test (dataProvider = "badGroupedFiles", expectedExceptions= SAMException.class)
    public void testBadGroupedFile(final String samFilename, final String fastq, final String secondEndFastq,
                                   final String [] groupFiles) throws IOException {
        final File samFile = new File(TEST_DATA_DIR,samFilename);
        final String tmpDir = IOUtil.getDefaultTmpDir().getAbsolutePath() + "/";
        final String [] args = new String[]{
                "--input", samFile.getAbsolutePath(),
                "--OUTPUT_PER_RG",
                "--OUTPUT_DIR", tmpDir,
        };
        runCommandLine(args);

        File f1;
        File f2;
        String fname1;
        String fname2;
        for(final String groupPUName : groupFiles)
        {
            fname1 = tmpDir + groupPUName + "_1.fastq";
            fname2 = tmpDir + groupPUName + "_2.fastq";
            f1 = new File(fname1);
            f2 = new File(fname2);
            f1.deleteOnExit();
            f1.deleteOnExit();
        }
    }

    @Test(dataProvider = "trimmedData")
    public void testTrimming(final String samFilename, final int read1Trim,
                             final int read1MaxBases, final int expectedRead1Length, final int read2Trim,
                             final int read2MaxBases, final int expectedRead2Length) throws IOException {

        final File samFile = new File(TEST_DATA_DIR, samFilename);
        final File pair1File = newTempFastqFile("pair1");
        final File pair2File = newTempFastqFile("pair2");

        convertFile(new String[]{
                "--input", samFile.getAbsolutePath(),
                "--FASTQ", pair1File.getAbsolutePath(),
                "--SECOND_END_FASTQ", pair2File.getAbsolutePath(),
                "--READ1_TRIM", Integer.toString(read1Trim),
                "--READ1_MAX_BASES_TO_WRITE", Integer.toString(read1MaxBases),
                "--READ2_TRIM", Integer.toString(read2Trim),
                "--READ2_MAX_BASES_TO_WRITE", Integer.toString(read2MaxBases)
        });

        for (final FastqRecord first : new FastqReader(pair1File)) {
            Assert.assertEquals(first.getReadString().length(), expectedRead1Length, "Incorrect read length");
            Assert.assertEquals(first.getBaseQualityString().length(), expectedRead1Length, "Incorrect quality string length");
        }
        for (final FastqRecord second : new FastqReader(pair2File)) {
            Assert.assertEquals(second.getReadString().length(), expectedRead2Length, "Incorrect read length");
            Assert.assertEquals(second.getBaseQualityString().length(), expectedRead2Length, "Incorrect quality string length");
        }
    }

    @DataProvider(name = "trimmedData")
    public Object[][] trimmedData() {
        return new Object[][] {
                // There are 13 bases in each of these reads
                {"ok/sorted-pair.sam", 6, 7, 7, 5, 8, 8}, // exact matches for everything
                {"ok/sorted-pair.sam", 7, 7, 6, 3, 8, 8}  // fewer or more bases
        };
    }

    private Set<String> createFastqReadHeaderSet(final File file) {
        final Set<String> set = new HashSet<>();
        final FastqReader freader = new FastqReader(file);
        while (freader.hasNext()) {
            final FastqRecord frec = freader.next();
            set.add(frec.getReadHeader());
        }
        return set ;
    }

    private Map<String,MatePair> createSamMatePairsMap(final File samFile) throws IOException {
        IOUtil.assertFileIsReadable(samFile);
        final Map<String, MatePair> map = new LinkedHashMap<>();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(samFile)) {

            for (final SAMRecord record : reader) {
                MatePair mpair = map.get(record.getReadName());
                if (mpair == null) {
                    mpair = new MatePair();
                    map.put(record.getReadName(), mpair);
                }
                mpair.add(record);
            }
        }
        return map;
    }


    private Map<String, Map<String, MatePair>> createPUPairsMap(final File samFile) throws IOException {
        IOUtil.assertFileIsReadable(samFile);
        final Map<String, Map<String, MatePair>> map = new LinkedHashMap<>();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(samFile)) {

            Map<String, MatePair> curFileMap;
            for (final SAMRecord record : reader) {
                final String platformUnit = record.getReadGroup().getPlatformUnit();
                curFileMap = map.get(platformUnit);
                if (curFileMap == null) {
                    curFileMap = new LinkedHashMap<>();
                    map.put(platformUnit, curFileMap);
                }

                MatePair mpair = curFileMap.get(record.getReadName());
                if (mpair == null) {
                    mpair = new MatePair();
                    curFileMap.put(record.getReadName(), mpair);
                }
                mpair.add(record);
            }
        }
        return map;
    }

    class MatePair {
        SAMRecord mate1 ;
        SAMRecord mate2 ;
        void add(final SAMRecord record) {
            if (!record.getReadPairedFlag()) throw new UserException("Record "+record.getReadName()+" is not paired");
            if (record.getFirstOfPairFlag()) {
                if (mate1 != null) throw new UserException("Mate 1 already set for record: "+record.getReadName());
                mate1 = record ;
            }
            else if (record.getSecondOfPairFlag()) {
                if (mate2 != null) throw new UserException("Mate 2 already set for record: "+record.getReadName());
                mate2 = record ;
            }
            else throw new UserException("Neither FirstOfPairFlag or SecondOfPairFlag is set for a paired record");
        }
    }

    private File newTempFastqFile(final String filename) throws IOException {
        if(filename == null) return null;
        final File file = BaseTest.createTempFile(filename, ".fastq");
        return file;
    }
}
