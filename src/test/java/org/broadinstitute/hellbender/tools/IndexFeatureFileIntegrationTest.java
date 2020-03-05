package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.linear.LinearIndex;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gtf.EnsemblGtfCodec;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfFeature;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public final class IndexFeatureFileIntegrationTest extends CommandLineProgramTest {

    final File ENSEMBL_GTF_TEST_FILE = new File(largeFileTestDir + "funcotator/ecoli_ds/gencode/ASM584v2/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.44.gtf");

    @Test
    public void testVCFIndex() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf");
        final File outName = createTempFile("test_variants_for_index.vcf", ".idx");

        final String[] args = {
                "-I", ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };

        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
        checkIndex(index, Arrays.asList("1", "2", "3", "4"));
    }

    @Test(groups={"bucket"})
    public void testVCFIndexOnCloud() throws IOException {
        final File testFile = getTestFile("test_variants_for_index.vcf");
        final String vcfOnGCS = BucketUtils.getTempFilePath(
                getGCPTestStaging() +"testIndexOnCloud", ".vcf");
        BucketUtils.copyFile(testFile.getAbsolutePath(), vcfOnGCS);

        final String[] args = new String[] {
                "IndexFeatureFile", "-I", vcfOnGCS
        };

        new Main().instanceMain(args);

        Assert.assertTrue(BucketUtils.fileExists(vcfOnGCS + ".idx"));

        final Index index = IndexFactory.loadIndex(vcfOnGCS + ".idx");
        Assert.assertTrue(index instanceof LinearIndex);
        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
        checkIndex(index, Arrays.asList("1", "2", "3", "4"));
    }

    @Test
    public void testVCFIndex_inferredName() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf");

        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
        };
        final Object res = this.runCommandLine(args);
        final Path tribbleIndex = Tribble.indexPath(ORIG_FILE.toPath());
        Assert.assertEquals(res, tribbleIndex.toAbsolutePath().toString());
        tribbleIndex.toFile().deleteOnExit();

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
        checkIndex(index, Arrays.asList("1", "2", "3", "4"));
    }

    @Test(expectedExceptions = UserException.NoSuitableCodecs.class)
    public void testIndexNonFeatureFileGZ() {
        final File ORIG_FILE = getTestFile("test_nonFeature_file.txt.blockgz.gz"); //made by bgzip
        final File outName = createTempFile("test_nonFeature_file.txt.blockgz.gz.", ".tbi");

        final String[] args = {
                "-I", ORIG_FILE.getAbsolutePath(),
                "-O", outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.NoSuitableCodecs.class)
    public void testIndexBCFFileGZ() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.bcf.blockgz.gz");  //made by bgzip
        final File outName = createTempFile("test_variants_for_index.bcf.blockgz.gz.", ".tbi");

        final String[] args = {
                "-I", ORIG_FILE.getAbsolutePath(),
                "-O", outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.class)
    public void testVCFGZIndex_tabixRequires_tbi_name() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf.blockgz.gz"); //made by bgzip
        final File outName = createTempFile("test_variants_for_index.blockgz.gz.", ".idx");

        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        this.runCommandLine(args);
    }

    @Test
    public void testVCFGZIndex_tabix() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf.blockgz.gz"); //made by bgzip
        final File outName = createTempFile("test_variants_for_index.blockgz.gz.",
            FileExtensions.TABIX_INDEX);

        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
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
    public void testVCFGZLargeHeaderIndex_tabix() throws IOException {
        // copy the input file, and create an index
        final File inputVCF = getTestFile("4featuresHG38Header.vcf.gz");
        final File tempDir = createTempDir("testVCFGZLargeHeaderIndex");
        final File inputCopy = new File(tempDir, inputVCF.getName());
        Files.copy(inputVCF.toPath(), inputCopy.toPath());
        final File outIndexFile = new File(tempDir, inputCopy.getName() + FileExtensions.TABIX_INDEX);

        final String[] args = {
                "-I" ,  inputCopy.getAbsolutePath(),
                "-O" ,  outIndexFile.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outIndexFile.getAbsolutePath());

        // use the location of every variant in the input as a query interval for the indexed copy of the same file
        try (final FeatureDataSource<VariantContext> originalSource = new FeatureDataSource<>(inputVCF);
                final FeatureDataSource<VariantContext> indexedSource = new FeatureDataSource<>(inputCopy)) {
            Iterator<VariantContext> originalIterator = originalSource.iterator();
            while (originalIterator.hasNext()) {
                final VariantContext originalVC = originalIterator.next();
                final Iterator<VariantContext> indexedIterator = indexedSource.query(new SimpleInterval(originalVC));
                Assert.assertTrue(indexedIterator.hasNext());
                final VariantContext queriedVC = indexedIterator.next();
                Assert.assertEquals(queriedVC.getContig(), originalVC.getContig());
                Assert.assertEquals(queriedVC.getStart(), originalVC.getStart());
                Assert.assertEquals(queriedVC.getEnd(), originalVC.getEnd());
                Assert.assertFalse(indexedIterator.hasNext());
            }
        }
    }

    @Test
    public void testVCFGZIndex_inferredName(){
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf.blockgz.gz"); //made by bgzip
        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
        };
        final Object res = this.runCommandLine(args);
        final File tabixIndex = new File(ORIG_FILE.getAbsolutePath() + FileExtensions.TABIX_INDEX);
        Assert.assertEquals(res, tabixIndex.getAbsolutePath());
        tabixIndex.deleteOnExit();

        Assert.assertTrue(tabixIndex.exists(), tabixIndex + " does not exists");
        final Index index = IndexFactory.loadIndex(tabixIndex.toString());
        Assert.assertTrue(index instanceof TabixIndex);

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
        checkIndex(index, Arrays.asList("1", "2", "3", "4"));
    }

    @Test(expectedExceptions = UserException.CouldNotIndexFile.class)
    public void testVCFGZIPIndex() throws IOException {
        //This tests blows up because the input file is not blocked gzipped
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf.gzip.gz"); //made by gzip
        final File outName = createTempFile("test_variants_for_index.gzip.gz.", ".tbi");
        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.CouldNotIndexFile.class)
    public void testVCFGZIPIndex_inferredName() throws IOException {
        //This tests blows up because the input file is not blocked gzipped
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf.gzip.gz"); //made by gzip
        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
        };
        final Object res = this.runCommandLine(args);
    }

    @Test
    public void testBCFIndex() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.bcf");
        final File outName = createTempFile("test_variants_for_index.bcf.", ".idx");

        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1"));
        checkIndex(index, Arrays.asList("1"));
    }

    // test disabled until https://github.com/samtools/htsjdk/issues/1323 is resolved
    @Test(enabled = false)
    public void testUncompressedBCF2_2Index() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.BCF22uncompressed.bcf");
        final File outName = createTempFile("test_variants_for_index.BCF22uncompressed.bcf", ".idx");

        final String[] args = {
                "-I", ORIG_FILE.getAbsolutePath(),
                "-O", outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.NoSuitableCodecs.class)
    public void testCompressedBCF2_2Index() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.BCF22compressed.bcf.blockgz.gz"); //made by bgzip
        final File outName = createTempFile("test_variants_for_index.BCF22compressed.bcf.blockgz.gz", ".idx");

        final String[] args = {
                "-I", ORIG_FILE.getAbsolutePath(),
                "-O", outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test
    public void testGVCFTreatedAsVCFIndex() {
        // Here we're testing what happens when we have a GVCF that is treated by the tool as a
        // regular VCF due to the lack of a .g.vcf extension
        final File ORIG_FILE = getTestFile("test_variants_for_index.gvcf_treated_as_vcf.vcf");
        final File outName = createTempFile("test_variants_for_index.gvcf_treated_as_vcf.vcf.", ".idx");

        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
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
        final File ORIG_FILE = getTestFile("test_variants_for_index.g.vcf");
        final File outName = createTempFile("test_variants_for_index.g.vcf.", ".idx");

        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
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

    private void testBedIndex(final File ORIG_FILE, final Class<? extends Index> indexClass) {
        final File outName = createTempFile(ORIG_FILE.getName(), (indexClass == TabixIndex.class) ?
            FileExtensions.TABIX_INDEX : ".idx");
        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(indexClass.isInstance(index));
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

    @Test
    public void testBedIndex() {
        testBedIndex(getTestFile("test_bed_for_index.bed"), LinearIndex.class);
    }

    @Test
    public void testBedGZIndex() {
        // made with bgzip
        testBedIndex(getTestFile("test_bed_for_index.bed.gz"), TabixIndex.class);
    }

    @Test
    public void testSAMPileupGZIndex() {
        final File ORIG_FILE = getTestFile("test_sampileup_for_index.pileup.gz"); // made with bgzip
        final File outName = createTempFile(ORIG_FILE.getName(), FileExtensions.TABIX_INDEX);

        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof TabixIndex);
        for (int chr: new int[]{1,2,3,4}) {  //note: unusual loop
            Assert.assertTrue(index.containsChromosome(String.valueOf(chr)), String.valueOf(chr));
        }
        for (int chr: new int[]{5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}) { //note: unusual loop
            Assert.assertFalse(index.containsChromosome(String.valueOf(chr)), String.valueOf(chr));
        }
        for (final String chr : Arrays.asList("X", "Y", "MT")){
            Assert.assertFalse(index.containsChromosome(chr), String.valueOf(chr));
        }

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testVCFIndex_missingFile() {
        final File ORIG_FILE = getTestFile("missing_file.vcf");
        final File outName = createTempFile("test_variants_for_index.vcf.", ".idx");

        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.CouldNotCreateOutputFile.class)
    public void testVCFIndex_cannotWrite() {
        final File ORIG_FILE = getTestFile("test_variants_for_index.vcf");
        final File tempDir = createTempDir("fred");
        final File doesNotExist = new File(tempDir, "bozo/");

        final File outName = new File(doesNotExist, "joe.txt");  //we can't write to this because parent does not exist

        final String[] args = {
                "-I" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    // Make sure we can index a VCF with 0 records without crashing:
    @Test
    public void testVCFWithNoRecords() {
        final File emptyVCF = getTestFile("header_only.vcf");
        final File output = createTempFile("header_only.vcf", ".idx");

        final String[] args = {
                "-I", emptyVCF.getAbsolutePath(),
                "-O", output.getAbsolutePath()
        };
        runCommandLine(args);

        Assert.assertTrue(output.exists());
        Assert.assertTrue(output.length() > 0);
    }

    @Test
    public void testEnsemblGtfIndex() {
        final File outName = createTempFile("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.44.gtf.", ".idx");

        final String[] args = {
                "-I" ,  ENSEMBL_GTF_TEST_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        Assert.assertEquals(index.getSequenceNames(), Collections.singletonList("Chromosome"));
        checkIndex(index, Collections.singletonList("Chromosome"));
    }

    @DataProvider
    Object[][] provideForTestEnsemblGtfIndexQuery() {
        return new Object[][] {
                {
                        new SimpleInterval("Chromosome", 3019160, 3020500),
                        1,
                        new SimpleInterval[] {new SimpleInterval("Chromosome", 3019161, 3020489)},
                        new String[] {"b2879"},
                        new String[] {"ssnA"},
                },
                {
                        new SimpleInterval("Chromosome", 3286269, 3288786),
                        4,
                        new SimpleInterval[] {
                                new SimpleInterval("Chromosome", 3285478, 3286269),
                                new SimpleInterval("Chromosome", 3286270, 3287025),
                                new SimpleInterval("Chromosome", 3287426, 3288010),
                                new SimpleInterval("Chromosome", 3288090, 3288785),
                        },
                        new String[] {
                                "b3140",
                                "b3141",
                                "b3142",
                                "b3143",
                        },
                        new String[] {
                                "agaD",
                                "agaI",
                                "yraH",
                                "yraI",
                        },
                }
        };
    }

    @Test(dataProvider = "provideForTestEnsemblGtfIndexQuery")
    public void testEnsemblGtfIndexQuery(final SimpleInterval interval,
                                         final Integer expectedNumResults,
                                         final SimpleInterval[] expectedFeatureIntervals,
                                         final String[] expectedGeneIds,
                                         final String[] expectedGeneNames) {
        // Test that we can query the file:
        try (final FeatureDataSource<GencodeGtfFeature> featureReader = new FeatureDataSource<>(ENSEMBL_GTF_TEST_FILE)) {
            final List<GencodeGtfFeature> features = featureReader.queryAndPrefetch(interval);

            Assert.assertEquals(features.size(), expectedNumResults.intValue());

            for ( int i = 0; i < expectedNumResults; ++i ) {
                final GencodeGtfFeature feature = features.get(i);
                Assert.assertEquals(feature.getGtfSourceFileType(), EnsemblGtfCodec.GTF_FILE_TYPE_STRING);
                Assert.assertEquals(feature.getChromosomeName(), expectedFeatureIntervals[ i ].getContig());
                Assert.assertEquals(feature.getStart(), expectedFeatureIntervals[ i ].getStart());
                Assert.assertEquals(feature.getEnd(), expectedFeatureIntervals[ i ].getEnd());

                Assert.assertEquals(feature.getGeneId(), expectedGeneIds[ i ]);
                Assert.assertEquals(feature.getGeneName(), expectedGeneNames[ i ]);
            }
        }
    }
}
