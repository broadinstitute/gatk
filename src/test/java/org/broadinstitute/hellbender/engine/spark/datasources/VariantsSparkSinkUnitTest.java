package org.broadinstitute.hellbender.engine.spark.datasources;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.utils.VCFHeaderReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.util.*;
import java.util.zip.GZIPInputStream;

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
                {hg19_chr1_1M_dbSNP, ".vcf", true},
                {hg19_chr1_1M_dbSNP, ".vcf.bgz", true},
                {hg19_chr1_1M_dbSNP, ".vcf.bgz", false}, // disable writing tabix
                {hg19_chr1_1M_dbSNP, ".vcf.gz", true},
                {hg19_chr1_1M_dbSNP, ".vcf.gz", false}, // disable writing tabix
                {hg19_chr1_1M_dbSNP_modified, ".vcf", true},
        };
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void variantsSinkTest(String vcf, String outputFileExtension, boolean writeTabixIndex) throws IOException {
        final File outputFile = createTempFile(outputFileName, outputFileExtension);
        assertSingleShardedWritingWorks(vcf, outputFile.getAbsolutePath(), writeTabixIndex);
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void variantsSinkHDFSTest(String vcf, String outputFileExtension, boolean writeTabixIndex) throws IOException {
        final String outputHDFSPath = MiniClusterUtils.getTempPath(cluster, outputFileName, outputFileExtension).toString();
        Assert.assertTrue(IOUtils.isHadoopUrl(outputHDFSPath));
        assertSingleShardedWritingWorks(vcf, outputHDFSPath, writeTabixIndex);
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void testWritingToAnExistingFileHDFS(String vcf, String outputFileExtension, boolean writeTabixIndex) throws IOException {
        final Path outputPath = MiniClusterUtils.getTempPath(cluster, outputFileName, outputFileExtension);
        final FileSystem fs = outputPath.getFileSystem(new Configuration());
        Assert.assertTrue(fs.createNewFile(outputPath));
        Assert.assertTrue(fs.exists(outputPath));
        assertSingleShardedWritingWorks(vcf, outputPath.toString(), writeTabixIndex);
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void testWritingToFileURL(String vcf, String outputFileExtension, boolean writeTabixIndex) throws IOException {
        String outputUrl = "file://" + createTempFile(outputFileName, outputFileExtension).getAbsolutePath();
        assertSingleShardedWritingWorks(vcf, outputUrl, writeTabixIndex);
    }

    @DataProvider
    public static Object[][] brokenCases() {
        return new Object[][]{
                {false, ".bcf"},
                {false, ".bcf.gz"},
                {true, ".g.bcf"},
                {true, ".g.bcf.gz"}
        };
    }

    @Test(dataProvider = "brokenCases", expectedExceptions = UserException.UnimplementedFeature.class)
    public void testBrokenGVCFCasesAreDisallowed(boolean writeGvcf, String extension) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        VariantsSparkSink.writeVariants(ctx, createTempFile("test", extension).toString(), null,
                new VCFHeader(), writeGvcf, Arrays.asList(1, 2, 4, 5), 2, 1, false);
    }

    @DataProvider
    public Object[][] gvcfCases(){
        return new Object[][]{
                {false, ".vcf"},
                {false, ".vcf.gz"},
                {true, ".g.vcf"},
                {true, ".g.vcf.gz"}
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
        VariantsSparkSink.writeVariants(ctx, output.toString(), ctx.parallelize(vcs), getHeader(), writeGvcf, Arrays.asList(100), 2, 1, true);

        checkFileExtensionConsistentWithContents(output.toString(), true);

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

    private void assertSingleShardedWritingWorks(String vcf, String outputPath, boolean writeTabixIndex) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        JavaRDD<VariantContext> variants = variantsSparkSource.getParallelVariantContexts(vcf, null);
        if (variants.getNumPartitions() == 1) {
            variants = variants.repartition(3); // repartition to more than 1 partition
        }
        VCFHeader header = getHeader(vcf);

        VariantsSparkSink.writeVariants(ctx, outputPath, variants, header, writeTabixIndex);

        checkFileExtensionConsistentWithContents(outputPath, writeTabixIndex);

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
            actualVcf = BaseTest.createTempFile(vcf, ".gz");
            Files.copy(IOUtils.getPath(vcf), actualVcf.toPath());
        } else {
            actualVcf = new File(vcf);
        }
        return VariantsSparkSourceUnitTest.getSerialVariantContexts(actualVcf.getAbsolutePath());
    }

    private void checkFileExtensionConsistentWithContents(String outputPath, boolean writeTabixIndex) throws IOException {
        final String outputFile = IOUtils.makeFilePathAbsolute(outputPath);
        boolean blockCompressed = IOUtil.isBlockCompressed(IOUtils.getPath(outputFile));
        String vcfFormat = getVcfFormat(outputFile);
        if (outputFile.endsWith(".vcf")) {
            Assert.assertEquals("VCF", vcfFormat);
            Assert.assertFalse(blockCompressed);
        } else if (outputFile.endsWith(".vcf.gz") || outputFile.endsWith(".vcf.bgz")) {
            Assert.assertEquals("VCF", vcfFormat);
            Assert.assertTrue(blockCompressed);
            // check that a tabix index file is created or not
            boolean tabixExists = Files.exists(IOUtils.getPath(outputPath + TabixUtils.STANDARD_INDEX_EXTENSION));
            Assert.assertEquals(tabixExists, writeTabixIndex);
        } else if (outputFile.endsWith(".bcf")) {
            Assert.assertEquals("BCF", vcfFormat);
            Assert.assertFalse(blockCompressed);
        } else if (outputFile.endsWith(".bcf.gz")) {
            Assert.assertEquals("BCF", vcfFormat);
            Assert.assertTrue(blockCompressed);
        }
    }

    private static String getVcfFormat(String outputFile) throws IOException {
        try (InputStream in = Files.newInputStream(IOUtils.getPath(outputFile))) {
            BufferedInputStream bis = new BufferedInputStream(in); // so mark/reset is supported
            return inferFromUncompressedData(IOUtil.isGZIPInputStream(bis) ? new GZIPInputStream(bis) : bis);
        }
    }

    private static String inferFromUncompressedData(final InputStream in) throws IOException {
        final byte b = (byte)in.read();
        in.close();
        switch (b) {
            case 'B':  return "BCF";
            case '#':  return "VCF";
        }
        return null;
    }

}
