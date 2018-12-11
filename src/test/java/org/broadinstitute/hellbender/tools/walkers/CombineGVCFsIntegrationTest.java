package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.variant.utils.VCFHeaderReader;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.codec.digest.DigestUtils;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections.IteratorUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

public class CombineGVCFsIntegrationTest extends CommandLineProgramTest {
    private static final List<String> NO_EXTRA_ARGS = Collections.emptyList();
    private static final List<String> ATTRIBUTES_TO_IGNORE = Arrays.asList(
            "RAW_MQ"); //MQ data format and key have changed since GATK3
    private static final File NA12878_HG37 = new File(toolsTestDir + "haplotypecaller/expected.testGVCFMode.gatk4.g.vcf");


    private static <T> void assertForEachElementInLists(final List<T> actual, final List<T> expected, final BiConsumer<T, T> assertion) {
        Assert.assertEquals(actual.size(), expected.size(), "different number of elements in lists:\n"
                + actual.stream().map(Object::toString).collect(Collectors.joining("\n","actual:\n","\n"))
                +  expected.stream().map(Object::toString).collect(Collectors.joining("\n","expected:\n","\n")));
        for (int i = 0; i < actual.size(); i++) {

            assertion.accept(actual.get(i), expected.get(i));
        }
    }

    @DataProvider
    public Object[][] gvcfsToCombine() {
        return new Object[][]{
                // Simple Test, spanning deletions
                {new File[]{getTestFile("spanningDel.1.g.vcf"),getTestFile("spanningDel.2.g.vcf")}, getTestFile("spanningDeletionRestrictToStartExpected.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},
                // Simple Test, multiple spanning deletions for one file
                {new File[]{getTestFile("spanningDel.many.g.vcf")}, getTestFile("testMultipleSpanningDeletionsForOneSample.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},
                // Simple Test, spanning deletions for haploid data
                {new File[]{getTestFile("spanningDel.many.haploid.g.vcf")}, getTestFile("testMultipleSpanningDeletionsForOneSampleHaploid.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},
                // Simple Test, spanning deletions for tetraploid data
                {new File[]{getTestFile("spanningDel.many.tetraploid.g.vcf")}, getTestFile("testMultipleSpanningDeletionsForOneSampleTetraploid.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},
                //Test that spanning deletion alleles in mixed sites don't get extended
                {new File[]{getTestFile("spanningDeletionBaseExtensionTest1.vcf"), getTestFile("spanningDeletionBaseExtensionTest2.vcf")},
                        getTestFile("spanningDeletionBaseExtensionTestExpected.g.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},
                // Testing BasePairResolutionInputs
                {new File[]{getTestFile("gvcf.basepairResolution.vcf")}, getTestFile("testBasepairResolutionInput.vcf"), Arrays.asList("-A", "ClippingRankSumTest"), b37_reference_20_21},
                // Interval Test
                {new File[]{getTestFile("gvcfExample1.vcf"),getTestFile("gvcfExample2.vcf"),}, getTestFile("IntervalTest.vcf"), Arrays.asList(" -L ",  "20:69485-69791"), b37_reference_20_21},
                // convert-to-base-pair-resolution argument test
                {new File[]{getTestFile("gvcfExample1.vcf"),getTestFile("gvcfExample2.vcf"),}, getTestFile("convertToBasePairResolution.vcf"), Arrays.asList(" -L ",  "20:69485-69791", "--" + CombineGVCFs.BP_RES_LONG_NAME), b37_reference_20_21},
                // Testing the breakBands argument " -L 1:69485-69791 --break-bands-at-multiples-of 5"
                {new File[]{getTestFile("gvcfExample1.vcf"),getTestFile("gvcfExample2.vcf"),}, getTestFile("testBreakBandsArgumet.vcf"), Arrays.asList(" -L ",  "20:69485-69791", "--" + CombineGVCFs.BREAK_BANDS_LONG_NAME, "5", "-A", "ClippingRankSumTest"), b37_reference_20_21},
                // Testing mismatched reference bases
                {new File[]{getTestFile("combine-gvcf-wrong-ref-input1.vcf"),getTestFile("combine-gvcf-wrong-ref-input2.vcf"),}, getTestFile("testWrongReferenceBaseBugFix.vcf"), Arrays.asList("-A", "ClippingRankSumTest"), b37_reference_20_21},
                //Testing allele-specific annotations
                {new File[]{getTestFile("NA12878.AS.chr20snippet.g.vcf"), getTestFile("NA12892.AS.chr20snippet.g.vcf")}, getTestFile("testAlleleSpecificAnnotations.vcf"), Arrays.asList("-G", "Standard", "-G", "AS_Standard"), b37_reference_20_21},
                //Testing allele-specific annotations missing AS_Standard Group
                {new File[]{getTestFile("NA12878.AS.chr20snippet.g.vcf"), getTestFile("NA12892.AS.chr20snippet.g.vcf")}, getTestFile("testAlleleSpecificAnnotationsNoGroup.vcf"), Arrays.asList("-G", "Standard", "-G", "AS_Standard"), b37_reference_20_21},
                //Test that trailing reference blocks are emitted with correct intervals
                {new File[]{getTestFile("gvcfExample1WithTrailingReferenceBlocks.g.vcf"), getTestFile("gvcfExample2WithTrailingReferenceBlocks.g.vcf")}, getTestFile("gvcfWithTrailingReferenceBlocksExpected.g.vcf"), NO_EXTRA_ARGS, b38_reference_20_21},
                // same test as the previous one, except with a band multiple specified
                {new File[]{getTestFile("gvcfExample1WithTrailingReferenceBlocks.g.vcf"), getTestFile("gvcfExample2WithTrailingReferenceBlocks.g.vcf")},
                        getTestFile("gvcfWithTrailingReferenceBlocksBandedExpected.g.vcf"),
                        Arrays.asList("--" + CombineGVCFs.BREAK_BANDS_LONG_NAME, "2000000"),
                        b38_reference_20_21},
                {new File[]{NA12878_HG37, getTestFile("YRIoffspring.chr20snippet.g.vcf")}, getTestFile("newMQcalc.combined.g.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},
        };
    }


    /*
    This test is useful for testing changes in GATK4 versus different versions of GATK3.
    To use, set GATK3_PATH to point to a particular version of gatk, then enable this test and run.

    It will cache the gatk3 outputs in a folder called gatk3results, it does it's best to avoid reusing bad results by
    comparing the md5 of the gatk3 path, input file path, reference, and commandline, but it doesn't know about internal changes to files
    You have to manually delete the cache if you make changes inside the input files.

    The expected outputs based on gatk3's results will be put in a folder called expectedResults.
    These are overwritten during each test, the names are based on the name of the existing expected output file.

    This method should be removed after GenotypeGVCFs has been completely validated against GATK3.
     */
    @Test(dataProvider = "gvcfsToCombine", enabled = false)
    public void compareToGATK3(File[] inputs, File outputFile, List<String> extraArgs, String reference) throws IOException, NoSuchAlgorithmException {
        final String GATK3_PATH = "/Users/emeryj/hellbender/gsa-unstable/target/package/GenomeAnalysisTK.jar";
        final String params = GATK3_PATH + inputs[0].getAbsolutePath() + extraArgs.stream().collect(Collectors.joining()) + reference;
        final String md5 = DigestUtils.md5Hex(params);
        final File gatk3ResultsDir = new File("gatk3results");
        if(! gatk3ResultsDir.exists()){
            Assert.assertTrue(gatk3ResultsDir.mkdir());
        }
        final File gatk3Result = new File(gatk3ResultsDir, md5 + ".vcf");
        if (true) {
            List<String> gatk3Command = new ArrayList<>(
                    Arrays.asList("java", "-jar", GATK3_PATH, "-T", "CombineGVCFs"));
            for (File f: inputs) {
                gatk3Command.add("-V");
                gatk3Command.add(f.getAbsolutePath());
            }
            gatk3Command.add("-o");
            gatk3Command.add(gatk3Result.getAbsolutePath());
            gatk3Command.add("-R");
            gatk3Command.add(reference);
            gatk3Command.addAll(extraArgs);

            runProcess(new ProcessController(), gatk3Command.toArray(new String[gatk3Command.size()]));
        } else {
            System.out.println("Found precomputed gatk3Result");
        }

        assertVariantContextsMatch(Arrays.asList(inputs), gatk3Result, extraArgs, reference, ATTRIBUTES_TO_IGNORE);
    }

    @Test(dataProvider = "gvcfsToCombine")
    public void compareToGATK3ExpectedResults(File[] inputs, File outputFile, List<String> extraArgs, String reference) throws IOException, NoSuchAlgorithmException {
        assertVariantContextsMatch(Arrays.asList(inputs), outputFile, extraArgs, reference, ATTRIBUTES_TO_IGNORE);
    }

    public static void runProcess(ProcessController processController, String[] command) {
        final ProcessSettings prs = new ProcessSettings(command);
        prs.getStderrSettings().printStandard(true);
        prs.getStdoutSettings().printStandard(true);
        final ProcessOutput output = processController.exec(prs);
        Assert.assertEquals(output.getExitValue(), 0, "Process exited with non-zero value. Command: "+ Arrays.toString(command) + "\n");
    }


    public void assertVariantContextsMatch(List<File> inputs, File expected, List<String> extraArgs, String reference, List<String> attributesToIgnore) throws IOException {
        final VCFHeader header = getHeaderFromFile(expected);

        runCombineGVCFSandAssertSomething(inputs, expected, extraArgs, (a, e) -> {
            VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, attributesToIgnore, header);
        }, reference);
    }


    public void runCombineGVCFSandAssertSomething(List<File> inputs, File expected, List<String> additionalArguments, BiConsumer<VariantContext, VariantContext> assertion, String reference) throws IOException {
        final File output = createTempFile("combinegvcfs", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
                .addOutput(output);
        for (File input: inputs) {
            args.addArgument("V", input.getAbsolutePath());
        }

        // Handling a difference in syntax between GATK3 and GATK4 wrt. annotation groups
        additionalArguments = additionalArguments.stream().map(a -> a.contains("Standard") ? a + "Annotation" : a).collect(Collectors.toList());
        additionalArguments.forEach(args::add);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> expectedVC = getVariantContexts(expected);
        final List<VariantContext> actualVC = getVariantContexts(output);
        assertForEachElementInLists(actualVC, expectedVC, assertion);
    }

    /**
     * Returns a list of VariantContext records from a VCF file
     *
     * @param vcfFile VCF file
     * @return list of VariantContext records
     * @throws IOException if the file does not exist or can not be opened
     */
    @SuppressWarnings({"unchecked"})
    private static List<VariantContext> getVariantContexts(final File vcfFile) throws IOException {
        try(final FeatureDataSource<VariantContext> variantContextFeatureDataSource = new FeatureDataSource<>(vcfFile)) {
            return IteratorUtils.toList(variantContextFeatureDataSource.iterator());
        }
    }

    private static VCFHeader getHeaderFromFile(final File vcfFile) throws IOException {
        try (SeekablePathStream stream = new SeekablePathStream(vcfFile.toPath())) {
            VCFHeader header = VCFHeaderReader.readHeaderFrom(stream);
            return header;
        }
    }

    @Test
    public void testOneStartsBeforeTwoAndEndsAfterwards() throws Exception {
        final File output = createTempFile("combinegvcfs", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addOutput(output);
        args.addVCF(getTestFile("gvcfExample1.vcf"));
        args.addVCF(getTestFile("gvcfExample2.vcf"));
        args.add(" -L 20:69485-69509");

        runCommandLine(args);

        final List<VariantContext> allVCs = getVariantContexts(output);

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

    @Test()
    public void testTetraploidRun() throws IOException {
        final File output = createTempFile("combinegvcfs", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addOutput(output);
        args.addArgument("variant","sample1:"+getToolTestDataDir()+"tetraploid-gvcf-1.vcf");
        args.addArgument("variant","sample2:"+getToolTestDataDir()+"tetraploid-gvcf-2.vcf");
        args.addArgument("variant","sample3:"+getToolTestDataDir()+"tetraploid-gvcf-3.vcf");
        args.addArgument("intervals", getToolTestDataDir() + "tetraploid-gvcfs.intervals");

        runCommandLine(args);

        final List<VariantContext> expectedVC = getVariantContexts(getTestFile("tetraploidRun.GATK3.g.vcf"));
        final List<VariantContext> actualVC = getVariantContexts(output);
        final VCFHeader header = getHeaderFromFile(output);
        assertForEachElementInLists(actualVC, expectedVC, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, Arrays.asList(), header));

    }

    @Test
    public void testTwoSpansManyBlocksInOne() throws Exception {
        final File output = createTempFile("combinegvcfs", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addOutput(output);
        args.addVCF(getTestFile("gvcfExample1.vcf"));
        args.addVCF(getTestFile("gvcfExample2.vcf"));
        args.add(" -L 20:69512-69634");

        runCommandLine(args);

        final List<VariantContext> allVCs = getVariantContexts(output);

        Assert.assertEquals(allVCs.size(), 5);
    }

    // Ensuring that no exception is thrown and that the resulting VCF is empty
    @Test
    public void testNoDataInInterval() throws Exception {
        final File output = createTempFile("combinegvcfs", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addOutput(output);
        args.addVCF(getTestFile("gvcfExample1.vcf"));
        args.addVCF(getTestFile("gvcfExample2.vcf"));
        args.add(" -L 1:69512-69634");

        runCommandLine(args);

        final List<VariantContext> allVCs = getVariantContexts(output);

        Assert.assertEquals(allVCs.size(), 0);
    }


    @Test
    public void testOneHasAltAndTwoHasNothing() throws Exception {
        final File output = createTempFile("combinegvcfs", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addOutput(output);
        args.addVCF(getTestFile("gvcfExample1.vcf"));
        args.addVCF(getTestFile("gvcfExample2.vcf"));
        args.add(" -L 20:69511");

        runCommandLine(args);

        final List<VariantContext> allVCs = getVariantContexts(output);

        final VariantContext first = allVCs.get(0);
        Assert.assertEquals(first.getStart(), 69511);
        Assert.assertEquals(first.getEnd(), 69511);
        Assert.assertEquals(first.getGenotypes().size(), 2);
    }

    @Test
    public void testOneHasAltAndTwoHasRefBlock() throws Exception {
        final File output = createTempFile("combinegvcfs", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addOutput(output);
        args.addVCF(getTestFile("gvcfExample1.vcf"));
        args.addVCF(getTestFile("gvcfExample2.vcf"));
        args.add("--" + CombineGVCFs.IGNORE_VARIANTS_THAT_START_OUTSIDE_INTERVAL );
        args.add(" -L 20:69635");

        runCommandLine(args);

        final List<VariantContext> allVCs = getVariantContexts(output);

        Assert.assertEquals(allVCs.size(), 1);

        final VariantContext first = allVCs.get(0);
        Assert.assertEquals(first.getStart(), 69635);
        Assert.assertEquals(first.getEnd(), 69635);
        Assert.assertEquals(first.getNAlleles(), 3);
        Assert.assertEquals(first.getGenotypes().size(), 2);
    }

    @Test
    public void testOneHasDeletionAndTwoHasRefBlock() throws Exception {
        final File output = createTempFile("combinegvcfs", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addOutput(output);
        args.addVCF(getTestFile("gvcfExample1.vcf"));
        args.addVCF(getTestFile("gvcfExample2.vcf"));
        args.add("--" + CombineGVCFs.IGNORE_VARIANTS_THAT_START_OUTSIDE_INTERVAL );
        args.add(" -L 20:69772-69783");

        runCommandLine(args);

        final List<VariantContext> allVCs = getVariantContexts(output);

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
    // This test asserts that the behavior is rational in the case of whole genome gvcfs which have variants starting at
    // the first and ending on the last base of each contig.
    public void testStartChromosome() throws Exception {
        final File output = createTempFile("combinegvcfs", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addOutput(output);
        args.addVCF(getTestFile("gvcfExample1.fullchrom.2021.vcf"));
        args.addVCF(getTestFile("gvcfExample2.vcf"));

        runCommandLine(args);

        final List<VariantContext> allVCs = getVariantContexts(output);

        final VariantContext first = allVCs.get(0);
        Assert.assertEquals(first.getStart(), 1);
        Assert.assertEquals(first.getEnd(), 69490);
        Assert.assertEquals(first.getContig(), "20");
        Assert.assertEquals(first.getNAlleles(), 2);
        Assert.assertEquals(first.getGenotypes().size(), 2);

        final VariantContext last = allVCs.get(19);
        Assert.assertEquals(last.getStart(), 69792);
        Assert.assertEquals(last.getEnd(), 63025520);
        Assert.assertEquals(last.getContig(), "20");
        Assert.assertEquals(last.getNAlleles(), 2);
        Assert.assertEquals(last.getGenotypes().size(), 2);

        final VariantContext nextchrom = allVCs.get(20);
        Assert.assertEquals(nextchrom.getStart(), 1);
        Assert.assertEquals(nextchrom.getEnd(), 69490);
        Assert.assertEquals(nextchrom.getContig(), "21");
        Assert.assertEquals(nextchrom.getNAlleles(), 2);
        Assert.assertEquals(nextchrom.getGenotypes().size(), 2);

        final VariantContext nextcromlast = allVCs.get(allVCs.size()-1);
        Assert.assertEquals(nextcromlast.getStart(), 69784);
        Assert.assertEquals(nextcromlast.getEnd(), 48129895);
        Assert.assertEquals(nextcromlast.getContig(), "21");
        Assert.assertEquals(nextcromlast.getNAlleles(), 2);
        Assert.assertEquals(nextcromlast.getGenotypes().size(), 2);
    }

    /** Should throw a BadInput exception as combining GVCFs with MNPs is unsupported. */
    @Test(expectedExceptions={UserException.BadInput.class})
    public void testCombineGvcfsWithMnps() throws Exception {
        final File output = createTempFile("combinegvcfs", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21));
        args.addOutput(output);
        args.addVCF(getTestFile("mnp.g.vcf"));
        runCommandLine(args);
    }
}
