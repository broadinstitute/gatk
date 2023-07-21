package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.writers.ShardingVCFWriter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.List;

public class GatkToolIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

    @Test
    public void testSitesOnlyMode() {
        File out = createTempFile("GTStrippedOutput", "vcf");
        String[] args = new String[] {
                "-V",  TEST_DIRECTORY + "vcf_with_genotypes.vcf",
                "--" + StandardArgumentDefinitions.SITES_ONLY_LONG_NAME,
                "-O",
                out.getAbsolutePath()};
        runCommandLine(Arrays.asList(args), SelectVariants.class.getSimpleName());

        // Assert that the genotype field has been stripped from the file
        Pair<VCFHeader, List<VariantContext>> results = VariantContextTestUtils.readEntireVCFIntoMemory(out.getAbsolutePath());

        Assert.assertFalse(results.getLeft().hasGenotypingData());
        for (VariantContext v: results.getRight()) {
            Assert.assertFalse(v.hasGenotypes());
        }
    }

    @Test (expectedExceptions = java.lang.IllegalArgumentException.class)
    // test asserting that if the reference dictionary exists but is not valid we get a more helpful exception than a null pointer exception
    public void testBrokenReferenceDictionaryErrorMessage() throws IOException {
        File out = createTempFile("GTStrippedOutput", "vcf");

        Path refCopy = Files.copy(IOUtils.getPath(hg19_chr1_1M_Reference), createTempPath("reference", ".fasta"), StandardCopyOption.REPLACE_EXISTING);
        Path indexCopy = Files.copy(ReferenceSequenceFileFactory.getFastaIndexFileName(IOUtils.getPath(hg19_chr1_1M_Reference)), ReferenceSequenceFileFactory.getFastaIndexFileName(refCopy));
        File emptyDict = new File(ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(refCopy).toString());
        IOUtils.deleteOnExit(indexCopy);

        emptyDict.createNewFile();
        IOUtils.deleteOnExit(emptyDict.toPath());

        String[] args = new String[] {
                "-R", refCopy.toString(),
                "-I", TEST_DIRECTORY + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.bam",
                "-O", out.getAbsolutePath()
        };

        runCommandLine(Arrays.asList(args), Mutect2.class.getSimpleName());
    }

    @Test
    public void testSharding() {
        final String outDir = createTempDir("GTShardedOutput").getAbsolutePath();
        final String fileBase = "test";
        final String out = Paths.get(outDir, fileBase + FileExtensions.COMPRESSED_VCF).toString();
        final String[] args = new String[] {
                "-V",  TEST_DIRECTORY + "example_variants_withSequenceDict.vcf",
                "-R", hg19MiniReference,
                "--" + StandardArgumentDefinitions.MAX_VARIANTS_PER_SHARD_LONG_NAME, "10",
                "-O", out};
        runCommandLine(Arrays.asList(args), SelectVariants.class.getSimpleName());

        // 12 total records in the test input should create 2 vcf shards with 10 and 2 records
        final Path basePath = Paths.get(outDir, fileBase).toAbsolutePath();
        final String firstShard = ShardingVCFWriter.getShardFilename(basePath, 0);
        final String secondShard = ShardingVCFWriter.getShardFilename(basePath, 1);
        final Pair<VCFHeader, List<VariantContext>> firstResults = VariantContextTestUtils.readEntireVCFIntoMemory(firstShard);
        final Pair<VCFHeader, List<VariantContext>> secondResults = VariantContextTestUtils.readEntireVCFIntoMemory(secondShard);
        Assert.assertEquals(firstResults.getValue().size(), 10, "First shard has wrong number of records");
        Assert.assertEquals(secondResults.getValue().size(), 2, "Second shard has wrong number of records");
        Assert.assertTrue(Files.exists(Paths.get(firstShard + FileExtensions.COMPRESSED_VCF_INDEX)));
        Assert.assertTrue(Files.exists(Paths.get(secondShard + FileExtensions.COMPRESSED_VCF_INDEX)));
    }
}
