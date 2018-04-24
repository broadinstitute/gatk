package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.io.Files;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
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
import java.util.*;

import static org.testng.Assert.assertEquals;

public final class VariantsSparkSinkUnitTest extends GATKBaseTest {
    private static final String SAMPLE = "sample";
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
                {hg19_chr1_1M_dbSNP, ".vcf.gz"},
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


    @DataProvider
    public static Object[][] brokenGVCFCases() {
        return new Object[][]{
                {"g.vcf.gz"},
                {"g.bcf"},
                {"g.bcf.gz"}
        };
    }

    @Test(dataProvider = "brokenGVCFCases", expectedExceptions = UserException.UnimplementedFeature.class)
    public void testBrokenGVCFCasesAreDisallowed(String extension) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        VariantsSparkSink.writeVariants(ctx, createTempFile("test", extension).toString(), null,
                                        new VCFHeader(), true, Arrays.asList(1, 2, 4, 5), 2, 1 );
    }

    @DataProvider
    public Object[][] gvcfCases(){
        return new Object[][]{
                {true, ".g.vcf"},
                {false, ".vcf"},
                {false, ".vcf.gz"},
                {false, ".bcf"},
                {false, ".bcf.gz"},
                //  {true, "g.vcf.gz"},  TODO enable this when https://github.com/broadinstitute/gatk/issues/4274 is resolved
                //  {true, ".g.bcf"},    TODO enable these when https://github.com/broadinstitute/gatk/issues/4303 is resolved
                //  {true, ".g.bcf.gz"}
        };
    }

    @Test(dataProvider = "gvcfCases")
    public void testEnableDisableGVCFWriting(boolean writeGvcf, String extension) throws IOException {
        List<VariantContext> vcs = new ArrayList<>();
        for(int i = 1; i <= 10; i++) {
            final Allele A = Allele.create("A", true);
            final VariantContext vc = new VariantContextBuilder("hand crafted", "1", i, i,
                                                                Arrays.asList(A, Allele.NON_REF_ALLELE))
                    .genotypes(new GenotypeBuilder(SAMPLE).alleles(Arrays.asList(A, A)).DP(10).GQ(10).PL(new int[]{0, 60, 10}).make())
                    .make();
            vcs.add(vc);
        }

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final File output = createTempFile(outputFileName, extension);
        VariantsSparkSink.writeVariants(ctx, output.toString(), ctx.parallelize(vcs), getHeader(), writeGvcf, Arrays.asList(100), 2, 1 );

        new CommandLineProgramTest(){
            @Override
            public String getTestedToolName(){
                return IndexFeatureFile.class.getSimpleName();
            }
        }.runCommandLine(new String[]{"-F", output.getAbsolutePath()});

        final List<VariantContext> writtenVcs = readVariants(output.toString());
        //if we are actually writing a gvcf, all the variant blocks will be merged into a single homref block with
        Assert.assertEquals(writtenVcs.size(), writeGvcf ? 1 : 10);
        Assert.assertEquals(writtenVcs.stream().mapToInt(VariantContext::getStart).min().getAsInt(), 1);
        Assert.assertEquals(writtenVcs.stream().mapToInt(VariantContext::getEnd).max().getAsInt(), 10);

    }

    private static VCFHeader getHeader() {
        final Set<VCFHeaderLine> headerlines = new LinkedHashSet<>();
        VCFStandardHeaderLines.addStandardFormatLines(headerlines, true,
                                                      VCFConstants.GENOTYPE_KEY,
                                                      VCFConstants.GENOTYPE_QUALITY_KEY,
                                                      VCFConstants.GENOTYPE_PL_KEY, VCFConstants.DEPTH_KEY);
        final SAMSequenceDictionary dict = new SAMSequenceDictionary(
                Collections.singletonList(new SAMSequenceRecord("1", 100)));
        final VCFHeader header = new VCFHeader(headerlines, Collections.singleton(SAMPLE));
        header.setSequenceDictionary(dict);
        return header;
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
