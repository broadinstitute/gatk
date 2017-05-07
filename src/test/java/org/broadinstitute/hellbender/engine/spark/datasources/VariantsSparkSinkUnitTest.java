package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.io.Files;
import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

public final class VariantsSparkSinkUnitTest extends BaseTest {
    private final String outputFileName = getClass().getSimpleName();
    private MiniDFSCluster cluster;

    @BeforeClass(alwaysRun = true)
    private void setupMiniCluster() throws IOException {
        cluster = MiniClusterUtils.getMiniCluster();
    }

    @AfterClass(alwaysRun = true)
    private void shutdownMiniCluster() {
        MiniClusterUtils.stopCluster(cluster);
    }


    @DataProvider(name = "loadVariants")
    public Object[][] loadVariants() {
        return new Object[][]{
                {hg19_chr1_1M_dbSNP, ".vcf"},
                {hg19_chr1_1M_dbSNP, ".vcf.bgz"},
                {hg19_chr1_1M_dbSNP_modified, ".vcf"},
        };
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void variantsSinkTest(String vcf, String outputFileExtension) throws IOException {
        final File outputFile = createTempFile(outputFileName, outputFileExtension);
        assertSingleShardedWritingWorks(vcf, outputFile.getAbsolutePath());
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void variantsSinkHDFSTest(String vcf, String outputFileExtension) throws IOException {
        final String outputHDFSPath = MiniClusterUtils.getTempPath(cluster, outputFileName, outputFileExtension).toString();
        Assert.assertTrue(BucketUtils.isHadoopUrl(outputHDFSPath));
        assertSingleShardedWritingWorks(vcf, outputHDFSPath);
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void testWritingToAnExistingFileHDFS(String vcf, String outputFileExtension) throws IOException {
        final Path outputPath = MiniClusterUtils.getTempPath(cluster, outputFileName, outputFileExtension);
        final FileSystem fs = outputPath.getFileSystem(new Configuration());
        Assert.assertTrue(fs.createNewFile(outputPath));
        Assert.assertTrue(fs.exists(outputPath));
        assertSingleShardedWritingWorks(vcf, outputPath.toString());
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void testWritingToFileURL(String vcf, String outputFileExtension) throws IOException {
        String outputUrl = "file://" + createTempFile(outputFileName, outputFileExtension).getAbsolutePath();
        assertSingleShardedWritingWorks(vcf, outputUrl);
    }

    private void assertSingleShardedWritingWorks(String vcf, String outputPath) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        JavaRDD<VariantContext> variants = variantsSparkSource.getParallelVariantContexts(vcf, null);
        if (variants.getNumPartitions() == 1) {
            variants = variants.repartition(3); // repartition to more than 1 partition
        }
        VCFHeader header = getHeader(vcf);

        VariantsSparkSink.writeVariants(ctx, outputPath, variants, header);

        JavaRDD<VariantContext> variants2 = variantsSparkSource.getParallelVariantContexts(outputPath, null);
        final List<VariantContext> writtenVariants = variants2.collect();

        VariantContextTestUtils.assertEqualVariants(readVariants(vcf), writtenVariants);
    }

    private VCFHeader getHeader(String vcf) throws IOException {
        final java.nio.file.Path vcfPath = IOUtils.getPath(vcf);
        try (SeekableStream stream = new SeekablePathStream(vcfPath)) {
            return VCFHeaderReader.readHeaderFrom(stream);
        } catch (IOException e) {
            throw new UserException("Failed to read VCF header from " + vcf + "\n Caused by:" + e.getMessage(), e);
        }
    }

    private List<VariantContext> readVariants(String vcf) throws IOException {
        File actualVcf;
        // work around TribbleIndexedFeatureReader not reading header from .bgz files
        if (vcf.endsWith(".bgz")) {
            actualVcf = File.createTempFile(vcf, ".gz");
            actualVcf.deleteOnExit();
            Files.copy(new File(vcf), actualVcf);
        } else {
            actualVcf = new File(vcf);
        }
        return VariantsSparkSourceUnitTest.getSerialVariantContexts(actualVcf.getAbsolutePath());
    }
}
