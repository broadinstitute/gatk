package org.broadinstitute.hellbender.tools.walkers;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import htsjdk.variant.variantcontext.VariantContext;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class CombineGVCFsIntegrationTest extends WalkerTest {
    public static String baseTestString(String args) {
        return "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -V "
                + privateTestDir + "gvcfExample1.vcf -V " + privateTestDir + "gvcfExample2.vcf" + args;
    }

    @Test
    public void testOneStartsBeforeTwoAndEndsAfterwards() throws Exception {
        final String cmd = baseTestString(" -L 1:69485-69509");
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
        spec.disableShadowBCF();
        final File gVCF = executeTest("testOneStartsBeforeTwoAndEndsAfterwards", spec).first.get(0);
        final List<VariantContext> allVCs = GATKVCFUtils.readVCF(gVCF).getSecond();

        Assert.assertEquals(allVCs.size(), 2, "Observed: " + allVCs);

        final VariantContext first = allVCs.get(0);
        Assert.assertEquals(first.getStart(), 69491);
        Assert.assertEquals(first.getEnd(), 69497);
        Assert.assertEquals(first.getGenotypes().size(), 2);
        Assert.assertTrue(first.getGenotype("NA1").isNoCall());
        Assert.assertTrue(first.getGenotype("NA2").isNoCall());

        final VariantContext second = allVCs.get(1);
        Assert.assertEquals(second.getStart(), 69498);
        Assert.assertEquals(second.getEnd(), 69506);
        Assert.assertEquals(second.getGenotypes().size(), 2);
        Assert.assertTrue(second.getGenotype("NA1").isNoCall());
        Assert.assertTrue(second.getGenotype("NA2").isNoCall());
    }

    @Test(enabled = true)
    public void testTetraploidRun() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -V:sample1 " + privateTestDir + "tetraploid-gvcf-1.vcf" +
                        " -V:sample2 " + privateTestDir + "tetraploid-gvcf-2.vcf" +
                        " -V:sample3 " + privateTestDir + "tetraploid-gvcf-3.vcf" +
                        " -L " + privateTestDir + "tetraploid-gvcfs.intervals",
                1,
                Arrays.asList("8472fcd43a41e7501b4709f0f1f5f432"));
        executeTest("combineSingleSamplePipelineGVCF", spec);
    }

    @Test(enabled= true)
    public void testMixedPloidyRun() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -V:sample1 " + privateTestDir + "haploid-gvcf-1.vcf" +
                        " -V:sample2 " + privateTestDir + "tetraploid-gvcf-2.vcf" +
                        " -V:sample3 " + privateTestDir + "diploid-gvcf-3.vcf" +
                        " -L " + privateTestDir + "tetraploid-gvcfs.intervals",
                1,
                Arrays.asList("f88b2b4d285276130cf088e7a03ca6a7"));
        executeTest("combineSingleSamplePipelineGVCF", spec);
    }

    @Test
    public void testTwoSpansManyBlocksInOne() throws Exception {
        final String cmd = baseTestString(" -L 1:69512-69634");
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
        spec.disableShadowBCF();
        final File gVCF = executeTest("testTwoSpansManyBlocksInOne", spec).first.get(0);
        final List<VariantContext> allVCs = GATKVCFUtils.readVCF(gVCF).getSecond();

        Assert.assertEquals(allVCs.size(), 5);
    }

    @Test
    public void testOneHasAltAndTwoHasNothing() throws Exception {
        final String cmd = baseTestString(" -L 1:69511");
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
        spec.disableShadowBCF();
        final File gVCF = executeTest("testOneHasAltAndTwoHasNothing", spec).first.get(0);
        final List<VariantContext> allVCs = GATKVCFUtils.readVCF(gVCF).getSecond();

        Assert.assertEquals(allVCs.size(), 1);

        final VariantContext first = allVCs.get(0);
        Assert.assertEquals(first.getStart(), 69511);
        Assert.assertEquals(first.getEnd(), 69511);
        Assert.assertEquals(first.getGenotypes().size(), 2);
    }

    @Test
    public void testOneHasAltAndTwoHasRefBlock() throws Exception {
        final String cmd = baseTestString(" -L 1:69635");
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
        spec.disableShadowBCF();
        final File gVCF = executeTest("testOneHasAltAndTwoHasRefBlock", spec).first.get(0);
        final List<VariantContext> allVCs = GATKVCFUtils.readVCF(gVCF).getSecond();

        Assert.assertEquals(allVCs.size(), 1);

        final VariantContext first = allVCs.get(0);
        Assert.assertEquals(first.getStart(), 69635);
        Assert.assertEquals(first.getEnd(), 69635);
        Assert.assertEquals(first.getNAlleles(), 3);
        Assert.assertEquals(first.getGenotypes().size(), 2);
    }

    @Test
    public void testOneHasDeletionAndTwoHasRefBlock() throws Exception {
        final String cmd = baseTestString(" -L 1:69772-69783");
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
        spec.disableShadowBCF();
        final File gVCF = executeTest("testOneHasDeletionAndTwoHasRefBlock", spec).first.get(0);
        final List<VariantContext> allVCs = GATKVCFUtils.readVCF(gVCF).getSecond();

        Assert.assertEquals(allVCs.size(), 3);

        final VariantContext first = allVCs.get(0);
        Assert.assertEquals(first.getStart(), 69772);
        Assert.assertEquals(first.getEnd(), 69776);
        Assert.assertEquals(first.getNAlleles(), 3);
        Assert.assertEquals(first.getGenotypes().size(), 2);

        final VariantContext second = allVCs.get(1);
        Assert.assertEquals(second.getStart(), 69773);
        Assert.assertEquals(second.getEnd(), 69774);
        Assert.assertEquals(second.getGenotypes().size(), 2);

        final VariantContext third = allVCs.get(2);
        Assert.assertEquals(third.getStart(), 69775);
        Assert.assertEquals(third.getEnd(), 69783);
        Assert.assertEquals(third.getGenotypes().size(), 2);
    }

    @Test
    public void testMD5s() throws Exception {
        final String cmd = baseTestString(" -L 1:69485-69791");
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("e1a888e8116cf59d53ad919634a18e6c"));
        spec.disableShadowBCF();
        executeTest("testMD5s", spec);
    }

    @Test
    public void testBasepairResolutionOutput() throws Exception {
        final String cmd = baseTestString(" -L 1:69485-69791 --convertToBasePairResolution");
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("e7e86722a49ad9730743c4952cdbedc7"));
        spec.disableShadowBCF();
        executeTest("testBasepairResolutionOutput", spec);
    }

    @Test
    public void testBreakBlocks() throws Exception {
        final String cmd = baseTestString(" -L 1:69485-69791 --breakBandsAtMultiplesOf 5");
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("bd279625ccfb3adcd39d97f07f3a236e"));
        spec.disableShadowBCF();
        executeTest("testBreakBlocks", spec);
    }

    @Test
    public void testSpanningDeletions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "spanningDel.1.g.vcf -V " + privateTestDir + "spanningDel.2.g.vcf",
                1,
                Arrays.asList("b22238e1ff584a157335429309fbfc5b"));
        spec.disableShadowBCF();
        executeTest("testSpanningDeletions", spec);
    }

    @Test
    public void testMultipleSpanningDeletionsForOneSample() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "spanningDel.many.g.vcf",
                1,
                Arrays.asList("b828ba5b69422ce32e234ea24d2df4c7"));
        spec.disableShadowBCF();
        executeTest("testMultipleSpanningDeletionsForOneSample", spec);
    }

    @Test
    public void testMultipleSpanningDeletionsForOneSampleHaploid() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "spanningDel.many.haploid.g.vcf",
                1,
                Arrays.asList("e707335ebd61bbe20775f76ad9b8c20d"));
        spec.disableShadowBCF();
        executeTest("testMultipleSpanningDeletionsForOneSampleHaploid", spec);
    }

    @Test
    public void testMultipleSpanningDeletionsForOneSampleTetraploid() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
                        " -V " + privateTestDir + "spanningDel.many.tetraploid.g.vcf",
                1,
                Arrays.asList("d4c22bd32d136414bfd7a6ebc5152026"));
        spec.disableShadowBCF();
        executeTest("testMultipleSpanningDeletionsForOneSampleTetraploid", spec);
    }

    @Test
    public void testWrongReferenceBaseBugFix() throws Exception {
        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -V " + (privateTestDir + "combine-gvcf-wrong-ref-input1.vcf"
                + " -V " + (privateTestDir + "combine-gvcf-wrong-ref-input2.vcf") + " -o %s --no_cmdline_in_header");
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("5fec22c9c8a0063f43c86ac86bb12e27"));
        spec.disableShadowBCF();
        executeTest("testWrongReferenceBaseBugFix",spec);

    }

    @Test
    public void testBasepairResolutionInput() throws Exception {
        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -V " + privateTestDir + "gvcf.basepairResolution.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("a839f81014758d1bdc900d59d35dd5bc"));
        spec.disableShadowBCF();
        executeTest("testBasepairResolutionInput", spec);
    }

    @Test
    public void testAlleleSpecificAnnotations() throws Exception {
        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -G Standard -G AS_Standard -V "
                + privateTestDir + "NA12878.AS.chr20snippet.g.vcf -V " + privateTestDir + "NA12892.AS.chr20snippet.g.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("26606fbfc3e9e5a8813307227d898b9a"));
        spec.disableShadowBCF();
        executeTest("testAlleleSpecificAnnotations", spec);
    }

    @Test
    public void testMissingAlleleSpecificAnnotationGroup() throws IOException {
        final File logFile = createTempFile("testMissingAlleleSpecificAnnotationGroup.log", ".tmp");
        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -V "
                + privateTestDir + "NA12878.AS.chr20snippet.g.vcf -V " + privateTestDir + "NA12892.AS.chr20snippet.g.vcf -log " +
                logFile.getAbsolutePath();
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList(""));
        spec.disableShadowBCF();
        executeTest("testMissingAlleleSpecificAnnotationGroup", spec);
        Assert.assertTrue(FileUtils.readFileToString(logFile).contains(ReferenceConfidenceVariantContextMerger.ADD_AS_STANDARD_MSG));
    }

    @Test
    public void testASMateRankSumAnnotation() throws Exception {
        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -G Standard -G AS_Standard -A AS_MQMateRankSumTest -V "
                + privateTestDir + "NA12878.AS.MateRankSum.chr20snippet.g.vcf -V " + privateTestDir + "NA12892.AS.MateRankSum.chr20snippet.g.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("9fd72e0f1f0f08a5aff6875225eec567"));
        spec.disableShadowBCF();
        executeTest("testASMateRankSumAnnotation", spec);
    }

    @Test
    public void testASInsertSizeRankSumAnnotation() throws Exception {
        final String cmd = "-T CombineGVCFs -R " + b37KGReference + " -o %s --no_cmdline_in_header -G Standard -G AS_Standard -A AS_InsertSizeRankSum -V "
                + privateTestDir + "NA12878.AS.InsertSizeRankSum.chr20snippet.g.vcf -V " + privateTestDir + "NA12892.AS.InsertSizeRankSum.chr20snippet.g.vcf";
        final WalkerTestSpec spec = new WalkerTestSpec(cmd, 1, Arrays.asList("b8e10684010702d11ee71538b1f201ac"));
        spec.disableShadowBCF();
        executeTest("testASInsertSizeRankSumAnnotation", spec);
    }


}
