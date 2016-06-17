package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.HashSet;
import java.util.Set;

public class HaplotypeCallerIntegrationTest extends CommandLineProgramTest {

    public static final String TEST_FILES_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/haplotypecaller/";

    /*
     * Test that in VCF mode we're consistent with past GATK4 results
     */
    @Test
    public void testVCFModeIsConsistentWithPastResults() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testVCFMode.gatk4.vcf");

        final String[] args = new String[] {
            "-I", NA12878_20_21_WGS_bam,
            "-R", b37_reference_20_21,
            "-L", "20:10000000-10100000",
            "-O", output.getAbsolutePath(),
            "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    /*
     * Test that in VCF mode we're >= 99% concordant with GATK3.4 results
     */
    @Test
    public void testVCFModeIsConcordantWithGATK3Results() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConcordantWithGATK3Results", ".vcf");
        final File gatk3Output = new File(TEST_FILES_DIR, "expected.testVCFMode.gatk3.4.vcf");

        final String[] args = new String[] {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.4 in VCF mode is < 99%");
    }

    /*
     * Test that in GVCF mode we're consistent with past GATK4 results
     */
    @Test
    public void testGVCFModeIsConsistentWithPastResults() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConsistentWithPastResults", ".g.vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testGVCFMode.gatk4.g.vcf");

        final String[] args = new String[] {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    /*
     * Test that in GVCF mode we're >= 99% concordant with GATK3 results
     */
    @Test
    public void testGVCFModeIsConcordantWithGATK3Results() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConcordantWithGATK3Results", ".g.vcf");
        final File gatk3Output = new File(TEST_FILES_DIR, "expected.testGVCFMode.gatk3.4.g.vcf");

        final String[] args = new String[] {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.4 in GVCF mode is < 99%");
    }

    @Test
    public void testBamoutProducesReasonablySizedOutput() {
        Utils.resetRandomGenerator();

        // We will test that when running with -bamout over the testInterval, we produce
        // a bam with a number of reads that is within 10% of what GATK3.4 produces with
        // -bamout over the same interval. This is just to test that we produce a reasonably-sized
        // bam for the region, not to validate the haplotypes, etc. We don't want
        // this test to fail unless there is a likely problem with -bamout itself (eg., empty
        // or truncated bam).
        final String testInterval = "20:10000000-10010000";
        final int gatk3BamoutNumReads = 5170;

        final File vcfOutput = createTempFile("testBamoutProducesReasonablySizedOutput", ".vcf");
        final File bamOutput = createTempFile("testBamoutProducesReasonablySizedOutput", ".bam");

        final String[] args = new String[] {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", testInterval,
                "-O", vcfOutput.getAbsolutePath(),
                "-bamout", bamOutput.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        try ( final ReadsDataSource bamOutReadsSource = new ReadsDataSource(bamOutput) ) {
            int actualBamoutNumReads = 0;
            for ( final GATKRead read : bamOutReadsSource ) {
                ++actualBamoutNumReads;
            }

            final int readCountDifference = Math.abs(actualBamoutNumReads - gatk3BamoutNumReads);
            Assert.assertTrue(((double)readCountDifference / gatk3BamoutNumReads) < 0.10,
                               "-bamout produced a bam with over 10% fewer/more reads than expected");
        }
    }

    /*
     * Calculate rough concordance between two vcfs, comparing only the positions, alleles, and the first genotype.
     */
    private double calculateConcordance( final File actual, final File expected ) {
        final Set<String> actualVCFKeys = new HashSet<>();
        final Set<String> expectedVCFKeys = new HashSet<>();
        int concordant = 0;
        int discordant = 0;

        try ( final FeatureDataSource<VariantContext> actualSource = new FeatureDataSource<>(actual, new VCFCodec());
              final FeatureDataSource<VariantContext> expectedSource = new FeatureDataSource<>(expected, new VCFCodec()) ) {

            for ( final VariantContext vc : actualSource ) {
                actualVCFKeys.add(keyForVariant(vc));
            }

            for ( final VariantContext vc : expectedSource ) {
                expectedVCFKeys.add(keyForVariant(vc));
            }

            for ( final String vcKey : actualVCFKeys ) {
                if ( ! expectedVCFKeys.contains(vcKey) ) {
                    ++discordant;
                }
                else {
                    ++concordant;
                }
            }

            for ( final String vcKey : expectedVCFKeys ) {
                if ( ! actualVCFKeys.contains(vcKey) ) {
                    ++discordant;
                }
            }
        }

        return (double)concordant / (double)(concordant + discordant);
    }

    private String keyForVariant( final VariantContext variant ) {
        return String.format("%s:%d-%d %s %s", variant.getContig(), variant.getStart(), variant.getEnd(),
                variant.getAlleles(), variant.getGenotype(0).getGenotypeString(false));
    }
}
