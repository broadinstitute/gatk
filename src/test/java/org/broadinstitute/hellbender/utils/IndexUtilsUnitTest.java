package org.broadinstitute.hellbender.utils;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Set;
import java.util.stream.Collectors;

public final class IndexUtilsUnitTest extends GATKBaseTest {

    @DataProvider(name= "okFeatureFiles")
    public Object[][] okFeatureFiles() {
        return new Object[][] {
                { new File(getToolTestDataDir(), "test_variants_for_index.vcf")},
                { new File(getToolTestDataDir(), "test_variants_for_index.g.vcf")},
                { new File(getToolTestDataDir(), "test_bed_for_index.bed")},
        };
    }

    @Test(dataProvider = "okFeatureFiles")
    public void testLoadIndex(final File featureFile) throws Exception {
        final Index index = IndexUtils.loadTribbleIndex(featureFile);
        Assert.assertNotNull(index);
    }

    @DataProvider(name= "okFeatureFilesTabix")
    public Object[][] okFeatureFilesTabix() {
        return new Object[][] {
                { new File(getToolTestDataDir(), "test_variants_for_index.vcf.bgz")},
        };
    }

    @Test(dataProvider = "okFeatureFilesTabix")
    public void testLoadTabixIndex(final File featureFile) throws Exception {
        final Index index = IndexUtils.loadTabixIndex(featureFile);
        Assert.assertNotNull(index);
    }

    @DataProvider(name= "failTabixIndexFiles")
    public Object[][] failTabixIndexFiles() {
        return new Object[][] {
                { new File(getToolTestDataDir(), "test_variants_for_index.vcf")},
                { new File(getToolTestDataDir(), "test_variants_for_index.g.vcf")},
                { new File(getToolTestDataDir(), "test_bed_for_index.bed")},
        };
    }

    @Test(dataProvider = "failTabixIndexFiles")
    public void testFailLoadTabixIndex(final File featureFile) throws Exception {
        final Index index = IndexUtils.loadTabixIndex(featureFile);
        Assert.assertNull(index);
    }

    @DataProvider(name= "failTribbleIndexFiles")
    public Object[][] failTribbleIndexFiles() {
        return new Object[][] {
                { new File(getToolTestDataDir(), "test_variants_for_index.vcf.bgz")},
        };
    }

    @Test(dataProvider = "failTribbleIndexFiles")
    public void testFailLoadTribbleIndex(final File featureFile) throws Exception {
        final Index index = IndexUtils.loadTribbleIndex(featureFile);
        Assert.assertNull(index);
    }

    @Test
    public void testLoadIndexAcceptOldIndex() throws Exception {
        final File featureFile = new File(getToolTestDataDir(), "test_variants_for_index.newerThanIndex.vcf");
        final File featureFileIdx = new File(getToolTestDataDir(), "test_variants_for_index.newerThanIndex.vcf.idx");
        final File tmpDir = GATKBaseTest.createTempDir("testLoadIndexAcceptOldIndex");

        final Path tmpFeatureFilePath= tmpDir.toPath().resolve(featureFile.toPath().getFileName());
        final Path tmpFeatureFileIdxPath= tmpDir.toPath().resolve(featureFileIdx.toPath().getFileName());

        Files.copy(featureFile.toPath(), tmpFeatureFilePath);
        Files.copy(featureFileIdx.toPath(), tmpFeatureFileIdxPath);

        final File tmpVcf= tmpFeatureFilePath.toFile();
        Thread.sleep(1000L); //wait a second
        tmpFeatureFilePath.toFile().setLastModified(System.currentTimeMillis()); //touch the file but not the index
        final Index index = IndexUtils.loadTribbleIndex(tmpVcf);
        Assert.assertNotNull(index);
        //this should NOT blow up (files newer than indices are tolerated)
    }

    @Test
    public void testLoadIndex_noIndex() throws Exception {
        final File featureFile = new File(getToolTestDataDir(), "test_variants_for_index.noIndex.vcf");
        final Index index = IndexUtils.loadTribbleIndex(featureFile);
        Assert.assertNull(index);
    }

    @Test
    public void testCheckIndexModificationTime() throws Exception {
        final File vcf = new File(getToolTestDataDir(), "test_variants_for_index.vcf");
        final File vcfIdx = new File(getToolTestDataDir(), "test_variants_for_index.vcf.idx");
        final Index index = IndexFactory.loadIndex(vcfIdx.getAbsolutePath());
        IndexUtils.checkIndexVersionAndModificationTime(vcf, vcfIdx, index);//no blowup
    }

    @Test
    public void testCreateSequenceDictionaryFromTribbleIndex() throws Exception {
        final SAMSequenceDictionary dict = IndexUtils.createSequenceDictionaryFromFeatureIndex(new File(getToolTestDataDir(), "test_variants_for_index.vcf"));
        final Set<String> contigs = dict.getSequences().stream().map(s -> s.getSequenceName()).collect(Collectors.toSet());
        Assert.assertEquals(contigs, Sets.newHashSet("1", "2", "3", "4"));
    }

    @Test
    public void testCreateSequenceDictionaryFromTabixIndex() throws Exception {
        final SAMSequenceDictionary dict = IndexUtils.createSequenceDictionaryFromFeatureIndex(new File(getToolTestDataDir(), "test_variants_for_index.vcf.bgz"));
        final Set<String> contigs = dict.getSequences().stream().map(s -> s.getSequenceName()).collect(Collectors.toSet());
        Assert.assertEquals(contigs, Sets.newHashSet("1", "2", "3", "4"));
    }

    @Test
    public void testIsSequenceDictionaryFromIndexPositive() throws Exception {
        final SAMSequenceDictionary dict = IndexUtils.createSequenceDictionaryFromFeatureIndex(new File(getToolTestDataDir(), "test_variants_for_index.vcf"));
        Assert.assertTrue(IndexUtils.isSequenceDictionaryFromIndex(dict));
    }

    @Test
    public void testIsSequenceDictionaryFromIndexNegative() throws Exception {
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("1", 99));
        dict.addSequence(new SAMSequenceRecord("2", 99));
        Assert.assertFalse(IndexUtils.isSequenceDictionaryFromIndex(dict));
    }

}
