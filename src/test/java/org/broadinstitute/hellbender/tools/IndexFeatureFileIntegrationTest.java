package org.broadinstitute.hellbender.tools;

import htsjdk.tribble.Tribble;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.linear.LinearIndex;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.util.TabixUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public final class IndexFeatureFileIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testVCFIndex() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf");
        final File outName = BaseTest.createTempFile("test_variants_for_index.vcf", ".idx");

        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
        checkIndex(index, Arrays.asList("1", "2", "3", "4"));
    }


    @Test
    public void testVCFIndex_inferredName() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf");

        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
        };
        final Object res = this.runCommandLine(args);
        final File tribbleIndex = Tribble.indexFile(ORIG_FILE);
        Assert.assertEquals(res, tribbleIndex.getAbsolutePath());
        tribbleIndex.deleteOnExit();

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
        checkIndex(index, Arrays.asList("1", "2", "3", "4"));
    }

    @Test(expectedExceptions = UserException.NoSuitableCodecs.class)
    public void testIndexNonFeatureFileGZ() {
        final File ORIG_FILE = getTestFile("test_nonFeature_file.txt.blockgz.gz"); //made by bgzip
        final File outName = BaseTest.createTempFile("test_nonFeature_file.txt.blockgz.gz.", ".idx");

        final String[] args = {
                "--feature_file", ORIG_FILE.getAbsolutePath(),
                "-O", outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.NoSuitableCodecs.class)
    public void testIndexBCFFileGZ() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.bcf.blockgz.gz");  //made by bgzip
        final File outName = BaseTest.createTempFile("test_variants_for_index.bcf.blockgz.gz.", ".idx");

        final String[] args = {
                "--feature_file", ORIG_FILE.getAbsolutePath(),
                "-O", outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.class)
    public void testVCFGZIndex_tabixRequires_tbi_name() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf.blockgz.gz"); //made by bgzip
        final File outName = BaseTest.createTempFile("test_variants_for_index.blockgz.gz.", ".idx");

        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        this.runCommandLine(args);
    }

    @Test
    public void testVCFGZIndex_tabix() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf.blockgz.gz"); //made by bgzip
        final File outName = BaseTest.createTempFile("test_variants_for_index.blockgz.gz.", TabixUtils.STANDARD_INDEX_EXTENSION);

        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof TabixIndex);

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
        checkIndex(index, Arrays.asList("1", "2", "3", "4"));
    }

    @Test
    public void testVCFGZIndex_inferredName(){
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf.blockgz.gz"); //made by bgzip
        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
        };
        final Object res = this.runCommandLine(args);
        final File tabixIndex = new File(ORIG_FILE.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
        Assert.assertEquals(res, tabixIndex.getAbsolutePath());
        tabixIndex.deleteOnExit();

        Assert.assertTrue(tabixIndex.exists(), tabixIndex + " does not exists");
        final Index index = IndexFactory.loadIndex(tabixIndex.toString());
        Assert.assertTrue(index instanceof TabixIndex);

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
        checkIndex(index, Arrays.asList("1", "2", "3", "4"));
    }

    @Test(expectedExceptions = TribbleException.MalformedFeatureFile.class)
    public void testVCFGZIPIndex() throws IOException {
        //This tests blows up because the input file is not blocked gzipped
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf.gzip.gz"); //made by gzip
        final File outName = BaseTest.createTempFile("test_variants_for_index.gzip.gz.", ".idx");
        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test(expectedExceptions = TribbleException.MalformedFeatureFile.class)
    public void testVCFGZIPIndex_inferredName() throws IOException {
        //This tests blows up because the input file is not blocked gzipped
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf.gzip.gz"); //made by gzip
        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
        };
        final Object res = this.runCommandLine(args);
    }

    @Test
    public void testBCFIndex() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.bcf");
        final File outName = BaseTest.createTempFile("test_variants_for_index.bcf.", ".idx");

        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1"));
        checkIndex(index, Arrays.asList("1"));
    }

    @Test(expectedExceptions = TribbleException.InvalidHeader.class)
    public void testUncompressedBCF2_2Index() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.BCF22uncompressed.bcf");
        final File outName = BaseTest.createTempFile("test_variants_for_index.BCF22uncompressed.bcf", ".idx");

        final String[] args = {
                "--feature_file", ORIG_FILE.getAbsolutePath(),
                "-O", outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.NoSuitableCodecs.class)
    public void testCompressedBCF2_2Index() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.BCF22compressed.bcf.blockgz.gz"); //made by bgzip
        final File outName = BaseTest.createTempFile("test_variants_for_index.BCF22compressed.bcf.blockgz.gz", ".idx");

        final String[] args = {
                "--feature_file", ORIG_FILE.getAbsolutePath(),
                "-O", outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test
    public void testGVCF_VCFIndex() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.gvcf.vcf");
        final File outName = BaseTest.createTempFile("test_variants_for_index.gvcf.vcf.", ".idx");

        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1"));
        checkIndex(index, Arrays.asList("1"));
    }

    @Test
    public void testGVCFIndex() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.gvcf");
        final File outName = BaseTest.createTempFile("test_variants_for_index.gvcf.", ".idx");

        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1"));
        checkIndex(index, Arrays.asList("1"));
    }

    private void checkIndex(Index index, List<String> chroms) {
        for (final String chrom : chroms) {
            Assert.assertTrue(index.containsChromosome(chrom));
        }
    }

    @Test
    public void testBedIndex() {
        final File ORIG_FILE = getTestFile("test_bed_for_index.bed");
        final File outName = BaseTest.createTempFile("test_bed_for_index.bed.", ".idx");

        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        for (int chr: new int[]{1,2,4}) {  //note: unusual loop
            Assert.assertTrue(index.containsChromosome(String.valueOf(chr)), String.valueOf(chr));
        }
        for (int chr: new int[]{3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}) { //note: unusual loop
            Assert.assertFalse(index.containsChromosome(String.valueOf(chr)), String.valueOf(chr));
        }
        for (final String chr : Arrays.asList("X", "Y", "MT")){
            Assert.assertFalse(index.containsChromosome(chr), String.valueOf(chr));
        }

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "4"));
    }


    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testVCFIndex_missingFile() {
        final File ORIG_FILE = getTestFile("missing_file.vcf");
        final File outName = BaseTest.createTempFile("test_variants_for_index.vcf.", ".idx");

        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.CouldNotCreateOutputFile.class)
    public void testVCFIndex_cannotWrite() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf");
        final File tempDir = BaseTest.createTempDir("fred");
        final File doesNotExist = new File(tempDir, "bozo/");

        final File outName = new File(doesNotExist, "joe.txt");  //we can't write to this because parent does not exist

        final String[] args = {
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }
}
