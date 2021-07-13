package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by David Benjamin on 2/17/17.
 */
@Test(groups = {"variantcalling"})
public class CreateSomaticPanelOfNormalsIntegrationTest extends CommandLineProgramTest {

    private static final File PON_VCFS_DIR = new File(toolsTestDir, "mutect/createpon/");

    // positions 10,000,000 - 11,000,000 of chr 20 and with most annotations removed
    private static final File GNOMAD = new File(largeFileTestDir, "very-small-gnomad.vcf");

    void checkPonVariants(final List<VariantContext> ponVariants, final List<VariantContext> inputVariants) {
        for (int n = 0; n < ponVariants.size(); n++) {
            final VariantContext ponVC = ponVariants.get(n);

            final VariantContext inputVC = inputVariants.get(n);

            Assert.assertEquals(ponVC.getAlleles(), inputVC.getAlleles());

            // by construction, every pon sample has each variant
            final double ponFraction = ponVC.getAttributeAsDouble(CreateSomaticPanelOfNormals.FRACTION_INFO_FIELD, 0.0);
            Assert.assertEquals(ponFraction, 1.0, 1.0e-6);

            List<Double> ponBeta = ponVC.getAttributeAsDoubleList(CreateSomaticPanelOfNormals.BETA_SHAPE_INFO_FIELD, 1.0);
            final double ponAF = ponBeta.get(0) / (ponBeta.get(0) + ponBeta.get(1));

            final int[] AD = inputVC.getGenotype(0).getAD();
            final double empiricalAF = (double) (MathUtils.sum(AD) - AD[0]) / MathUtils.sum(AD);

            Assert.assertEquals(ponAF, empiricalAF, 0.2);
        }
    }

    @Test
    public void testTwoIdenticalVcfs() throws IOException {
        final File vcf1 = new File(PON_VCFS_DIR, "sample1.vcf");
        final File vcf2 = new File(PON_VCFS_DIR, "sample1-copy.vcf");

        final List<File> inputs = Arrays.asList(vcf1, vcf2);

        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(inputs, new SimpleInterval("20",1,63_025_520));
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        final File output = createTempFile("ponvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(new File(b37Reference))
                .add("V", genomicsDBUri)
                .addOutput(output);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> inputVariants = VariantContextTestUtils.streamVcf(inputs.get(0)).collect(Collectors.toList());
        final List<VariantContext> ponVariants = VariantContextTestUtils.streamVcf(output).collect(Collectors.toList());

        // this checks that even filtered variants make it into the pon
        Assert.assertEquals(ponVariants.size(), inputVariants.size());

        checkPonVariants(ponVariants, inputVariants);

        // Test again by requesting BCF codec stream from GenomicsDBFeatureReader
        args.add(GenomicsDBArgumentCollection.USE_BCF_CODEC_LONG_NAME, true);
        Utils.resetRandomGenerator();
        runCommandLine(args);
        final List<VariantContext> ponVariantsWithBCFCodecStreaming = VariantContextTestUtils.streamVcf(output).collect(Collectors.toList());
        Assert.assertEquals(ponVariantsWithBCFCodecStreaming.size(), inputVariants.size());
        checkPonVariants(ponVariantsWithBCFCodecStreaming, inputVariants);
    }

    @Test
    public void testJustOneVcf() throws IOException {
        final File vcf = new File(PON_VCFS_DIR, "sample1.vcf");

        final List<File> inputs = Arrays.asList(vcf);

        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(inputs, new SimpleInterval("20",1,63_025_520));
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        final File output = createTempFile("ponvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(new File(b37Reference))
                .add("V", genomicsDBUri)
                .addOutput(output);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> ponVariants = VariantContextTestUtils.streamVcf(output).collect(Collectors.toList());
        Assert.assertTrue(ponVariants.isEmpty());

        // Test again by requesting BCF codec stream from GenomicsDBFeatureReader
        args.add(GenomicsDBArgumentCollection.USE_BCF_CODEC_LONG_NAME, true);
        Utils.resetRandomGenerator();
        runCommandLine(args);
        Assert.assertTrue(VariantContextTestUtils.streamVcf(output).collect(Collectors.toList()).isEmpty());
    }


    @Test
    public void testGnomad() throws IOException {
        final File vcf = new File(PON_VCFS_DIR, "sample1.vcf");

        final List<File> inputs = Arrays.asList(vcf);

        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(inputs, new SimpleInterval("20",10_000_000,11_000_000));
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        final File output = createTempFile("ponvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(new File(b37Reference))
                .add("V", genomicsDBUri)
                .add(CreateSomaticPanelOfNormals.MIN_SAMPLE_COUNT_LONG_NAME, "1")
                .add(M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, GNOMAD.getAbsolutePath())
                .addOutput(output);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> ponVariants = VariantContextTestUtils.streamVcf(output).collect(Collectors.toList());
        Assert.assertTrue(!ponVariants.isEmpty());

        // here is a variant that should definitely be skipped as germline
        Assert.assertTrue(ponVariants.stream().mapToInt(VariantContext::getStart).noneMatch(start -> start == 10_000_117));

        // Test again by requesting BCF codec stream from GenomicsDBFeatureReader
        args.add(GenomicsDBArgumentCollection.USE_BCF_CODEC_LONG_NAME, true);
        Utils.resetRandomGenerator();
        runCommandLine(args);
        final List<VariantContext> ponVariantsFromBCFCodecStream = VariantContextTestUtils.streamVcf(output).collect(Collectors.toList());
        Assert.assertTrue(!ponVariantsFromBCFCodecStream.isEmpty());
        Assert.assertTrue(ponVariantsFromBCFCodecStream.stream().mapToInt(VariantContext::getStart).noneMatch(start -> start == 10_000_117));
    }

    void checkPonVariantsFromThreeVcfs(final List<VariantContext> ponVariants) {
        Assert.assertEquals(ponVariants.size(), 3);

        final double[] ponFractions = ponVariants.stream()
                .mapToDouble(vc -> vc.getAttributeAsDouble(CreateSomaticPanelOfNormals.FRACTION_INFO_FIELD, 0.0))
                .toArray();
        Assert.assertEquals(ponFractions[0], 2.0/3, 0.01);
        Assert.assertEquals(ponFractions[1], 1, 0.01);
        Assert.assertEquals(ponFractions[2], 1, 0.01);

        final List<BetaDistributionShape> betas = ponVariants.stream().map(vc -> {
            final List<Double> ponBeta = vc.getAttributeAsDoubleList(CreateSomaticPanelOfNormals.BETA_SHAPE_INFO_FIELD, 1.0);
            return new BetaDistributionShape(ponBeta.get(0), ponBeta.get(1));
        }).collect(Collectors.toList());
    }

    /**
     * Test on several samples over a few sites
     * 20:577114  G,A       ADs: 100,1  100,0   100,3
     * 20:577180  G,GAAAACA ADs: 10,30  20,20   30,10
     * 20:577548  C,G       ADs: 1000,10    10000,100   10000,100
     */
    @Test
    public void testThreeVcfs() throws IOException {
        final File vcf2 = new File(PON_VCFS_DIR, "sample2.vcf");
        final File vcf3 = new File(PON_VCFS_DIR, "sample3.vcf");
        final File vcf4 = new File(PON_VCFS_DIR, "sample4.vcf");

        final List<File> inputs = Arrays.asList(vcf2,vcf3,vcf4);

        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(inputs, new SimpleInterval("20",1,63_025_520));
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        final File output = createTempFile("ponvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(new File(b37Reference))
                .add("V", genomicsDBUri)
                .addOutput(output);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> ponVariants = VariantContextTestUtils.streamVcf(output).collect(Collectors.toList());
        checkPonVariantsFromThreeVcfs(ponVariants);

        // Test again by requesting BCF codec stream from GenomicsDBFeatureReader
        args.add(GenomicsDBArgumentCollection.USE_BCF_CODEC_LONG_NAME, true);
        Utils.resetRandomGenerator();
        runCommandLine(args);
        checkPonVariantsFromThreeVcfs(VariantContextTestUtils.streamVcf(output).collect(Collectors.toList()));
    }
}