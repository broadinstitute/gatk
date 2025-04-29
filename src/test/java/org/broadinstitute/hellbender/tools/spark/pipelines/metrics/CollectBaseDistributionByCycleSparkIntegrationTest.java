package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Integration tests for {@link CollectBaseDistributionByCycleSpark}.
 */
public final class CollectBaseDistributionByCycleSparkIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(
            "src/test/resources/org/broadinstitute/hellbender/metrics/analysis/CollectBaseDistributionByCycle");

    //Note: the 'expected' results in this test come from running picard 1.130
    //Note: these tests use the same data and results as the non-spark ones, by design
    //Note: we don't test the contents of the chart pdf

    @DataProvider(name = "CollectBaseDistributionByCycle")
    private Iterator<Object[]> makeCollectBaseDistributionByCycleData() {
        final List<Object[]> list = new ArrayList<>();

        list.add(new Object[]{"valid.bam",  null, "valid.CollectBaseDistributionByCycle.txt", true, false, false});
        list.add(new Object[]{"valid.cram", new File(TEST_DATA_DIR, "valid.fasta").getAbsolutePath(), "valid.CollectBaseDistributionByCycle.txt", true, false, false});

        list.add(new Object[]{"first5000a.bam", null, "CollectBaseDistributionByCycle.txt", true, false, false});
        list.add(new Object[]{"first5000a.cram", b37_reference_20_21, "CollectBaseDistributionByCycle.txt", true, false, false});

        list.add(new Object[]{"originalQuals.chr1.1-1K.bam", null, "CollectBaseDistributionByCycle.origQuals.txt", true, false, false});

        list.add(new Object[]{"example_pfFail_reads.bam", null, "CollectBaseDistributionByCycle.pfReads.txt", true, false, false});
        list.add(new Object[]{"example_pfFail_reads.bam", null, "CollectBaseDistributionByCycle.pfOnly.txt", true, true, false});

        list.add(new Object[]{"unmapped.bam", null, "CollectBaseDistributionByCycle.unmapped.ALIGNED_READS_ONLY_false.txt", true, false, false});
        list.add(new Object[]{"unmapped.bam", null, "CollectBaseDistributionByCycle.unmapped.ALIGNED_READS_ONLY_true.txt", false, false, true});

        return list.iterator();
    }

    @Test(dataProvider = "CollectBaseDistributionByCycle", groups = {"R", "spark"})
    public void test(final String unsortedBamName, final String referenceName, final String expectedFileName, final boolean makePdf, final boolean pfReadsOnly, final boolean alignedReadsOnly) throws IOException {
        final File unsortedBam = new File(TEST_DATA_DIR, unsortedBamName);
        final File reference = referenceName == null ? null : new File(referenceName);
        final File expectedFile = new File(TEST_DATA_DIR, expectedFileName);

        final File outfile = GATKBaseTest.createTempFile("test", ".metrics");
        final File pdf = GATKBaseTest.createTempFile("test", ".pdf");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.addRaw("-" + "I");
        args.addRaw(unsortedBam.getCanonicalPath());
        args.addRaw("-" + "O");
        args.addRaw(outfile.getCanonicalPath());
        if (makePdf) {
            args.addRaw("--" + "chart");
            args.addRaw(pdf.getCanonicalPath());
        }
        if (reference != null){
            args.addRaw("-" + "R");
            args.addRaw(reference.getCanonicalPath());
        }
        args.addRaw("--" + "pf-reads-only");
        args.addRaw(pfReadsOnly);
        args.addRaw("--" + "aligned-reads-only");
        args.addRaw(alignedReadsOnly);

        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }

    @Test
    public void testInGatkLiteDocker() throws IOException {
        final String gatkLiteDockerProperty = System.getProperty("IN_GATKLITE_DOCKER");

        try {
            System.setProperty("IN_GATKLITE_DOCKER", "true");

            final File unsortedBam = new File(TEST_DATA_DIR, "unsortedBamName");

            final File outfile = GATKBaseTest.createTempFile("test", ".metrics");
            final File pdf = GATKBaseTest.createTempFile("test", ".pdf");

            ArgumentsBuilder args = new ArgumentsBuilder();
            args.addRaw("-" + "I");
            args.addRaw(unsortedBam.getCanonicalPath());
            args.addRaw("-" + "O");
            args.addRaw(outfile.getCanonicalPath());
            args.addRaw("--" + "chart");
            args.addRaw(pdf.getCanonicalPath());
            args.addRaw("--" + "pf-reads-only");
            args.addRaw(false);
            args.addRaw("--" + "aligned-reads-only");
            args.addRaw(false);

            try {
                runCommandLine(args.getArgsArray());
                Assert.fail("Excepted RuntimeException for running in GATK Lite docker");
            }
            catch(final RuntimeException e) {
                Assert.assertTrue(e.getMessage().contains("Generating a chart file requires R, which is not available in the GATK Lite Docker image.")); 
            }

        }
        finally {
            if(gatkLiteDockerProperty != null) {
                System.setProperty("IN_GATKLITE_DOCKER", gatkLiteDockerProperty);
            }
            else{
                System.clearProperty("IN_GATKLITE_DOCKER");
            } 
        }

    }

    @Test
    public void testGetRScriptResource() {
        // Make sure the RScript resource can be resolved
        Assert.assertNotNull(CollectBaseDistributionByCycleSpark.getBaseDistributionByCycleRScriptResource());
    }
}
