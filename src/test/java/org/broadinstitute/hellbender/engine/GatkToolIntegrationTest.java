package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
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
        createTempFile("reference", ".dict"); // file because we want an empty file that "exists" for the reference
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
}
