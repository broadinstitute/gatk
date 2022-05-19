package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.CommandLineProgramTester;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class ReblockGVCFIntegrationTest extends CommandLineProgramTest {

    private static final String hg38_reference_20_21 = largeFileTestDir + "Homo_sapiens_assembly38.20.21.fasta";
    private static final String b37_reference_20_21 = largeFileTestDir + "human_g1k_v37.20.21.fasta";
    public static final String WARP_PROD_REBLOCKING_ARGS = " -do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40 ";

    @DataProvider(name = "getCommandLineArgsForExactTest")
    public Object[][] getCommandLineArgsForExactTest() {
        return new Object[][]{
                //covers inputs with old format "MQ" annotation
                {getTestFile("gvcfForReblocking.g.vcf"), getTestFile("testJustOneSample.expected.g.vcf"), " -L chr20:69771 -rgq-threshold 19", hg38_reference_20_21},
                //Broad production arguments on WGS data
                {getTestFile("prodWgsInput.g.vcf "), getTestFile("prodWgsOutput.g.vcf"), WARP_PROD_REBLOCKING_ARGS, hg38Reference},
                //Exome data with AS annotations and zero DP regression test
                {getTestFile("prodWesInput.g.vcf "), getTestFile("prodWesOutput.g.vcf"), WARP_PROD_REBLOCKING_ARGS, hg38Reference}
        };
    }

    @Test(dataProvider = "getCommandLineArgsForExactTest")
    public void testWithExactComparison(final File input, final File expected, final String extraArgs, final String reference) throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -O %s -R " + reference +
                        " -V " + input.getAbsolutePath() +
                        extraArgs +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(expected.getAbsolutePath()));
        spec.executeTest("testWithExactComparison", this);
    }

    @Test
    public void testGVCFReblockingIsContiguous() throws Exception {
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final File expected = new File(largeFileTestDir + "testProductionGVCF.expected.g.vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .add("V", largeFileTestDir + "NA12878.prod.chr20snippet.g.vcf.gz")
                .add("rgq-threshold", "20")
                .add("L", "20:60001-1000000")
                .add("A", "Coverage")
                .add("A", "RMSMappingQuality")
                .add("A", "ReadPosRankSumTest")
                .add("A", "MappingQualityRankSumTest")
                .add("disable-tool-default-annotations", true)
                .addOutput(output);
        runCommandLine(args);

        final CommandLineProgramTester validator = ValidateVariants.class::getSimpleName;
        final ArgumentsBuilder args2 = new ArgumentsBuilder();
        args2.add("R", b37_reference_20_21);
        args2.add("V", output.getAbsolutePath());
        args2.add("L", "20:60001-1000000");
        args2.addRaw("-gvcf");
        validator.runCommandLine(args2);  //will throw a UserException if GVCF isn't contiguous

        try (final FeatureDataSource<VariantContext> actualVcs = new FeatureDataSource<>(output);
             final FeatureDataSource<VariantContext> expectedVcs = new FeatureDataSource<>(expected)) {
            GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                    (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                            Collections.emptyList(), Collections.emptyList()));
        }
    }

    @Test
    public void testContiguityWithSpanningDels() {
        final File input = new File(getToolTestDataDir() + "regressionTests.vcf");
        final File intervals = new File(getToolTestDataDir() + "regressionContinuityIntervals.list");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .add("V", input.getAbsolutePath())
                .addOutput(output);
        runCommandLine(args);

        final CommandLineProgramTester validator = ValidateVariants.class::getSimpleName;
        final ArgumentsBuilder args2 = new ArgumentsBuilder();
        args2.add("R", hg38Reference);
        args2.add("V", output.getAbsolutePath());
        args2.add("L", intervals.getAbsolutePath());
        args2.addRaw("-gvcf");  //TODO: add --no-overlaps after ValidateVariants update is merged
        validator.runCommandLine(args2);  //will throw a UserException if GVCF isn't contiguous
    }

    @Test  //absolute minimal output
    public void testOneSampleAsForGnomAD() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-drop-low-quals -do-qual-approx -L chr20:69485-69791 -O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "gvcfForReblocking.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false" +
                        " -A Coverage -A RMSMappingQuality -A ReadPosRankSumTest -A MappingQualityRankSumTest --disable-tool-default-annotations true",
                Arrays.asList(getToolTestDataDir() + "testOneSampleAsForGnomAD.expected.g.vcf"));
        spec.executeTest("testOneSampleDropLows", this);
    }

    @Test  //covers non-ref AD and non-ref GT corrections
    public void testNonRefADCorrection() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "nonRefAD.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "testNonRefADCorrection.expected.g.vcf"));
        spec.executeTest("testNonRefADCorrection", this);
    }

    @Test //covers inputs with "RAW_MQ" annotation
    public void testRawMQInput() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "prod.chr20snippet.withRawMQ.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "prod.chr20snippet.withRawMQ.expected.g.vcf"));
        spec.executeTest("testRawMQInput", this);
    }

    @Test
    public void testASAnnotationsAndSubsetting() throws Exception {
        //some subsetting, but never dropping the first alt
        //also has multi-allelic that gets trimmed with ref block added
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + b37_reference_20_21 +
                        " -drop-low-quals -do-qual-approx -V " + "src/test/resources/org/broadinstitute/hellbender/tools/walkers/CombineGVCFs/NA12878.AS.chr20snippet.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "expected.NA12878.AS.chr20snippet.reblocked.g.vcf"));
        spec.executeTest("testASAnnotationsAndSubsetting", this);

        //one case where first alt is dropped
        final IntegrationTestSpec spec2 = new IntegrationTestSpec(
                "-O %s -R " + b37_reference_20_21 +
                        " -drop-low-quals -do-qual-approx -V " + "src/test/resources/org/broadinstitute/hellbender/tools/walkers/CombineGVCFs/NA12892.AS.chr20snippet.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "expected.NA12892.AS.chr20snippet.reblocked.g.vcf"));
        spec2.executeTest("testASAnnotationsAndSubsetting2", this);

        //big test for as we ran for gnomADv3
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", "src/test/resources/large/ReblockGVCF/spanDel.exome.chr20.vcf")
                .add("do-qual-approx", true)
                .add("drop-low-quals", true)
                .add("rgq-threshold", "10")
                .add("L", "chr20")
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        List<VariantContext> actual = VariantContextTestUtils.getVariantContexts(output);
        actual.stream().forEach(a -> {
            if (!a.getGenotype(0).isHomRef()) {
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_RAW_QUAL_APPROX_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_SB_TABLE_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(a, GATKVCFConstants.AS_VARIANT_DEPTH_KEY,
                    VCFHeaderLineCount.R);
        } });

    }

    @Test
    public void testNewCompressionScheme() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + b37_reference_20_21 +
                        " -drop-low-quals -do-qual-approx -V " + "src/test/resources/org/broadinstitute/hellbender/tools/walkers/CombineGVCFs/NA12878.AS.chr20snippet.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false" +
                        " --floor-blocks -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60",
                Arrays.asList(getToolTestDataDir() + "expected.NA12878.AS.chr20snippet.reblocked.hiRes.g.vcf"));
        spec.executeTest("testNewCompressionScheme", this);
    }

    @Test
    public void testAggressiveQualFiltering() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + hg38_reference_20_21 +
                        " -drop-low-quals -do-qual-approx -V " + getToolTestDataDir() + "gvcfForReblocking.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false" +
                        " --floor-blocks" +
                        " --" + GenotypeCalculationArgumentCollection.CALL_CONFIDENCE_LONG_NAME + " 65.0",
                Arrays.asList(getToolTestDataDir() + "expected.aggressiveQualFiltering.g.vcf"));
        spec.executeTest("testVariantQualFiltering", this);
    }

    @Test
    public void testMQHeadersAreUpdated() throws Exception {
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", getToolTestDataDir() + "justHeader.g.vcf")
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        Pair<VCFHeader, List<VariantContext>> actual = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        VCFHeader header = actual.getLeft();
        List<VCFInfoHeaderLine> infoLines = new ArrayList<>(header.getInfoHeaderLines());
        //check all the headers in case there's one old and one updated
        for (final VCFInfoHeaderLine line : infoLines) {
            if (line.getID().equals(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED)) {
                Assert.assertTrue(line.getType().equals(VCFHeaderLineType.Float));
                Assert.assertTrue(line.getDescription().contains("deprecated"));
            } else if (line.getID().equals(GATKVCFConstants.MAPPING_QUALITY_DEPTH_DEPRECATED)) {
                Assert.assertTrue(line.getDescription().contains("deprecated"));
            }
        }
    }

    @Test
    public void testReReblocking() {
        final File input = new File(getToolTestDataDir() + "alreadyReblocked.chr22snippet.vcf");
        final File output = createTempFile("rereblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", input)
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        final List<VariantContext> inputVCs = VariantContextTestUtils.readEntireVCFIntoMemory(input.getAbsolutePath()).getRight();
        final List<VariantContext> outputVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight();

        Assert.assertTrue(inputVCs.size() > outputVCs.size());
        Assert.assertTrue(outputVCs.size() == 19);
        //hom ref blocks change, but variants stay the same
        Assert.assertEquals(inputVCs.stream().filter(vc -> !vc.getGenotype(0).isHomRef()).count(), outputVCs.stream().filter(vc -> !vc.getGenotype(0).isHomRef()).count());
        List<String> inGenotypes= inputVCs.stream().filter(vc -> !vc.getGenotype(0).isHomRef()).map(vc -> vc.getGenotype(0)).map(Genotype::toString).collect(Collectors.toList());
        List<String> outGenotypes = outputVCs.stream().filter(vc -> !vc.getGenotype(0).isHomRef()).map(vc -> vc.getGenotype(0)).map(Genotype::toString).collect(Collectors.toList());
        Assert.assertTrue(inGenotypes.containsAll(outGenotypes)); //will check ref and alt alleles as part of genotype string representation
        Assert.assertTrue(outputVCs.get(18).isVariant());

        //all ref blocks have MIN_DP
        Assert.assertEquals(outputVCs.stream().filter(vc -> vc.getGenotype(0).hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY)).count(), outputVCs.size() - outGenotypes.size());
        //all variants have GQ
        Assert.assertEquals(outputVCs.stream().filter(vc -> vc.getGenotype(0).hasGQ()).count(), outputVCs.size());
        //we didn't ask to drop GQ0s, but they might get merged together
        Assert.assertEquals(inputVCs.stream().anyMatch(vc -> vc.getGenotype(0).getGQ() == 0), outputVCs.stream().anyMatch(vc -> vc.getGenotype(0).getGQ() == 0));
    }

    @Test
    public void testOverlappingDeletions() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + hg38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "overlappingDeletions.hc.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "expected.overlappingDeletions.g.vcf"));
        spec.executeTest("testOverlappingDeletions", this);

        final File output = createTempFile("reblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", getToolTestDataDir() + "overlappingDels2.vcf" )
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        //regression test for dropped span del allele
        Pair<VCFHeader, List<VariantContext>> actual = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final List<VariantContext> variants = actual.getRight();
        Assert.assertEquals(variants.size(), 3);
        Assert.assertEquals(variants.get(1).getGenotype(0).getAllele(0), Allele.SPAN_DEL);
    }

    @Test
    public void testHomRefCalls() throws IOException {
        final File input = new File(getToolTestDataDir() + "dropGQ0Dels.g.vcf");
        final File output = createTempFile("dropGQ0Dels.reblocked", ".g.vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", input)
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        final List<VariantContext> inputVCs = VariantContextTestUtils.readEntireVCFIntoMemory(input.getAbsolutePath()).getRight();
        final List<VariantContext> outputVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight();

        Assert.assertEquals(outputVCs.size(), 3);
        Assert.assertEquals(outputVCs.get(0).getStart(), inputVCs.get(0).getStart());
        Assert.assertEquals(outputVCs.get(outputVCs.size()-1).getEnd(), inputVCs.get(inputVCs.size()-1).getEnd());
        Assert.assertEquals(outputVCs.get(1).getGenotype(0).getGQ(), 0);  //there should be a GQ0 block in the middle from a crap variant in the input
    }

    @Test
    public void testMultipleInputs() {
        //run with multiple inputs split from chr20:19995000-19998999 of prod.chr20snippet.withRawMQ.g.vcf
        //note that an event is duplicated in shard1 and shard2 because it spans the boundary
        final File output = createTempFile("multi-input", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", getToolTestDataDir() + "chr20.shard3.g.vcf")
                .add("V", getToolTestDataDir() + "chr20.shard2.g.vcf")
                .add("V", getToolTestDataDir() + "chr20.shard1.g.vcf")
                .add("V", getToolTestDataDir() + "chr20.shard0.g.vcf")
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        final File output2 = createTempFile("single-input",".vcf");
        final ArgumentsBuilder args2 = new ArgumentsBuilder();
        args2.add("V", getToolTestDataDir() + "prod.chr20snippet.withRawMQ.g.vcf")
                .add("L", "chr20:19995000-19998999")
                .addReference(hg38Reference)
                .addOutput(output2);
        runCommandLine(args2);

        try (final FeatureDataSource<VariantContext> actualVcs = new FeatureDataSource<>(output);
             final FeatureDataSource<VariantContext> expectedVcs = new FeatureDataSource<>(output2)) {
            GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                    (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                            Collections.emptyList(), Collections.emptyList()));
        }
    }

    @Test(expectedExceptions = UserException.class)
    public void testMixedSamples() {
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", getToolTestDataDir() + "justHeader.g.vcf") //sample "Sample"
            .add("V", getToolTestDataDir() + "nonRefAD.g.vcf") //sample "HK017-0046"
            .addReference(hg38Reference)
            .addOutput(output);
        runCommandLine(args);
    }

    @Test
    //we had some external GVCFs that each went through CombineGVCFs for some reason, so GTs all went to ./.
    //chr4 is regression for overlapping ref blocks induced by a hom ref with a real alt
    //also test variants and ref blocks with no DP
    public void testNoCallGenotypes() {
        final File output = createTempFile("reblockedgvcf", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", getToolTestDataDir() + "noCallGTs.g.vcf")
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        Pair<VCFHeader, List<VariantContext>> actual = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        final List<VariantContext> variants = actual.getRight();
        final List<String> variantKeys = variants.stream().map(VariantContextTestUtils::keyForVariant).collect(Collectors.toList());
        final Map<String, VariantContext> resultMap = new LinkedHashMap<>();
        for (int i = 0; i < variants.size(); i++) {
            resultMap.put(variantKeys.get(i), variants.get(i));
        }

        final List<String> expectedHomVarKeys = Arrays.asList(
                "chr22:10514994-10514994 G*, [<NON_REF>, A]",
                "chr22:10515170-10515170 C*, [<NON_REF>, T]",
                "chr22:10515223-10515223 G*, [<NON_REF>, C]",
                "chrY:6067982-6067982 T*, [<NON_REF>, G]");

        final List<String> expectedHetKeys = Arrays.asList(
                "chr4:24339212-24339218 GTATATA*, [<NON_REF>, G]",
                "chr22:10515120-10515120 A*, [<NON_REF>, AAAGC]",
                "chr22:10515223-10515223 G*, [<NON_REF>, C]");

        final List<String> expectedHomRefKeys = Arrays.asList(
                "chr4:24339199-24339209 T*, [<NON_REF>]",
                "chr4:24339210-24339211 G*, [<NON_REF>]",
                "chr4:24339219-24339219 T*, [<NON_REF>]",
                "chr4:24339221-24339221 T*, [<NON_REF>]",
                "chr4:24339222-24339226 A*, [<NON_REF>]",
                "chr22:10515118-10515118 G*, [<NON_REF>]",
                "chrY:6067983-6068157 A*, [<NON_REF>]");

        Assert.assertTrue(variantKeys.containsAll(expectedHomVarKeys));
        Assert.assertTrue(variantKeys.containsAll(expectedHetKeys));
        Assert.assertTrue(variantKeys.containsAll(expectedHomRefKeys));
        Assert.assertTrue(variants.size() == 31);
    }

    @Test
    //here genotype call disagrees with PLs and subsetting can be wrong and lead to no confidence hom-ref if we rely on GT
    public void testPosteriorDisagreementNeedingSubset() {
        final File funkyDragenVariant = new File(getToolTestDataDir() + "HG002.snippet.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", funkyDragenVariant)
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        final Pair<VCFHeader, List<VariantContext>> outVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());
        Assert.assertEquals(outVCs.getRight().size(), 2);  //one variant and one ref block from trimmed deletion
        final VariantContext testResult = outVCs.getRight().get(0);
        final Genotype g = testResult.getGenotype(0);
        Assert.assertFalse(testResult.hasAttribute(VCFConstants.END_KEY));
        Assert.assertTrue(g.isHetNonRef());
        Assert.assertEquals(testResult.getAlternateAlleles().size(), 3);
        Assert.assertTrue(testResult.getAlternateAlleles().contains(Allele.NON_REF_ALLELE));
        Assert.assertEquals(testResult.getReference().getBaseString().length(), 1);  //alleles are properly trimmed

        final VariantContext newRefBlock = outVCs.getRight().get(1);
        Assert.assertTrue(newRefBlock.hasAttribute(VCFConstants.END_KEY));
        Assert.assertEquals(newRefBlock.getAttributeAsInt(VCFConstants.END_KEY, 0), 103392934);
        final Genotype refG = newRefBlock.getGenotype(0);
        Assert.assertTrue(refG.isHomRef());
        Assert.assertEquals(newRefBlock.getAlternateAlleles().size(), 1);
        Assert.assertTrue(newRefBlock.getAlternateAlleles().contains(Allele.NON_REF_ALLELE));
    }

    @Test
    public void testLeftContigBoundary() {
        final File funkyDragenVariant = new File(getToolTestDataDir() + "testLeftContigBoundary.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", funkyDragenVariant)
                .addReference(hg38Reference)
                .addOutput(output);
        runCommandLine(args);

        final List<VariantContext> outVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight();
        Assert.assertEquals(outVCs.size(), 2);
        Assert.assertEquals(outVCs.get(0).getStart(), 1);
        Assert.assertEquals(outVCs.get(1).getStart(), 1);
    }

    @Test
    public void testTreeScoreThreshold() {
        final File inputGvcf = new File(getToolTestDataDir() + "treeScoreGvcf.g.vcf");
        final File output = createTempFile("reblockedgvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", inputGvcf)
                .addOutput(output)
                .addReference(hg38_reference_20_21)
                .add(ReblockGVCF.TREE_SCORE_THRESHOLD_LONG_NAME, 0.3)
                .add(ReblockGVCF.ANNOTATIONS_TO_KEEP_LONG_NAME, "TREE_SCORE")
                .add("do-qual-approx", true);

        runCommandLine(args);

        final List<VariantContext> outVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight();
        for (VariantContext vc : outVCs) {
            if(vc.getStart() == 69511) {
                Assert.assertEquals(vc.getGenotype(0).getGQ(), 99, "Site chr20:69511 should not have been changed from GQ 99.");
            }
            if(vc.getStart() == 69512) {
                Assert.assertEquals(vc.getGenotype(0).getGQ(), 99, "Ref block chr20:69512 should not have been changed from GQ 99.");
            }
            if(vc.getStart() == 69767) {
                Assert.assertEquals(vc.getGenotype(0).getGQ(), 0, "Ref block chr20:69767 should have been grouped with GQ 0 block.");
                Assert.assertEquals(vc.getEnd(), 69783, "Ref block chr20:69767 should have been expanded to include site 69771.");
            }
            Assert.assertNotEquals(vc.getStart(), 69771, "Site chr20:69771 should have been made into a ref block due to low TREE_SCORE.");
            if(vc.getStart() == 69785){
                Assert.assertEquals(vc.getGenotype(0).getGQ(), 0, "Site chr20:69785 should have been made GQ 0 due to missing TREE_SCORE");
            }
            if(!vc.isReferenceBlock()) {
                Assert.assertTrue(vc.hasAttribute(GATKVCFConstants.TREE_SCORE));
            }
        }
    }
}
