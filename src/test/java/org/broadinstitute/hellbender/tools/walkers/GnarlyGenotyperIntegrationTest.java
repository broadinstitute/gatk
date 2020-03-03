package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.testutils.VariantContextTestUtils.getVariantContexts;

/**
 * Created by gauthier on 3/9/18.
 */
public class GnarlyGenotyperIntegrationTest extends CommandLineProgramTest {
    private static final List<String> NO_EXTRA_ARGS = Collections.emptyList();

    @DataProvider(name="VCFdata")
    public Object[][] getVCFdata() {
        return new Object[][]{
                // Simple Test, spanning deletions; standard calling confidence
                //No variants outside requested intervals; no SNPs with QUAL < 60, no INDELs with QUAL < 69?; has star alleles after deletion at chr20:263497; has AC, AF, AN, DP, ExcessHet, FS, MQ, (MQRankSum), (ReadPosRankSum), SOR, QD; has called genotypes
                {new File[]{getTestFile("sample1.vcf"), getTestFile("sample2.vcf"), getTestFile("sample3.vcf"), getTestFile("sample4.vcf"), getTestFile("sample5.vcf")},
                         getTestFile("fiveSampleTest.vcf"), null, Arrays.asList(new SimpleInterval("chr20", 251370, 252000), new SimpleInterval("chr20", 263000, 265600)), Arrays.asList("--merge-input-intervals", "--only-output-calls-starting-in-intervals"), b38_reference_20_21},
                // Same as above, but with GenomicsDB using VCFCodec for interchange
                {new File[]{getTestFile("sample1.vcf"), getTestFile("sample2.vcf"), getTestFile("sample3.vcf"), getTestFile("sample4.vcf"), getTestFile("sample5.vcf")},
                         getTestFile("fiveSampleTest.vcf"), null, Arrays.asList(new SimpleInterval("chr20", 251370, 252000), new SimpleInterval("chr20", 263000, 265600)), Arrays.asList("--merge-input-intervals", "--only-output-calls-starting-in-intervals", "--genomicsdb-use-vcf-codec"), b38_reference_20_21},


                //lower calling confidence
                //same as above except (different intervals and) with SNPs with 40 < QUAL < 60 and INDELs with 49 < QUAL < 69
                {new File[]{getTestFile("sample1.vcf"), getTestFile("sample2.vcf"), getTestFile("sample3.vcf"), getTestFile("sample4.vcf"), getTestFile("sample5.vcf")},
                         getTestFile("fiveSampleTest.lowerCallThreshold.vcf"), null, Arrays.asList(new SimpleInterval("chr20", 250865, 348163)), Arrays.asList("-stand-call-conf 10"), b38_reference_20_21},
                // Same as above, but with GenomicsDB using VCFCodec for interchange
                {new File[]{getTestFile("sample1.vcf"), getTestFile("sample2.vcf"), getTestFile("sample3.vcf"), getTestFile("sample4.vcf"), getTestFile("sample5.vcf")},
                         getTestFile("fiveSampleTest.lowerCallThreshold.vcf"), null, Arrays.asList(new SimpleInterval("chr20", 250865, 348163)), Arrays.asList("-stand-call-conf 10", "--genomicsdb-use-vcf-codec"), b38_reference_20_21},

                //using latest reblocking output with allele-specific annotations
                //has all of the above plus AS_AltDP, AS_FS, AS_MQ, AS_MQRankSum, AS_QD, AS_ReadPosRankSum
                {new File[]{new File(getToolTestDataDir() + "/../variantutils/ReblockGVCF/expected.NA12878.AS.chr20snippet.reblocked.g.vcf"),
                         new File(getToolTestDataDir() + "/../variantutils/ReblockGVCF/expected.NA12892.AS.chr20snippet.reblocked.g.vcf")},
                         getTestFile("twoSampleAS.vcf"), getTestFile("twoSampleASDB.vcf"), Arrays.asList(new SimpleInterval("20")), NO_EXTRA_ARGS, b37_reference_20_21},
                // Same as above, but with GenomicsDB using VCFCodec for interchange
                {new File[]{new File(getToolTestDataDir() + "/../variantutils/ReblockGVCF/expected.NA12878.AS.chr20snippet.reblocked.g.vcf"),
                        new File(getToolTestDataDir() + "/../variantutils/ReblockGVCF/expected.NA12892.AS.chr20snippet.reblocked.g.vcf")},
                        getTestFile("twoSampleAS.vcf"), getTestFile("twoSampleASDB.vcf"), Arrays.asList(new SimpleInterval("20")), Arrays.asList("--genomicsdb-use-vcf-codec"), b37_reference_20_21},

                //using legacy reblocking data with no raw GT count values
                //has ExcessHet, calculated from genotypes
                {new File[]{new File(getToolTestDataDir() + "noGTCount.sample1.chr20snippet.vcf"),
                        new File(getToolTestDataDir() + "noGTCount.sample2.chr20snippet.vcf")},
                        getTestFile("noGTCount.expected.chr20snippet.vcf"), null, Arrays.asList(new SimpleInterval("chr20")), NO_EXTRA_ARGS, hg38Reference
                }
        };
    }

    /**
     * Have a customized checker as the VariantContextUtils.assertEqualVariants version does a
     * simplistic ToString() comparison that does not account for differences in Scientific versus Decimal
     * notations.
     * @param actualVC List of actual variant contexts
     * @param expectedVC List of expected variant contexts
     */
    private static void assertEqualVariants(
            final List<VariantContext> actualVC, final List<VariantContext> expectedVC) {
        Utils.nonNull(actualVC, "actualVC");
        Utils.nonNull(expectedVC, "expectedVC");
        Assert.assertEquals(actualVC.size(), expectedVC.size());
        for (int i = 0; i < actualVC.size(); i++) {
            VariantContextTestUtils.assertVariantContextsAreEqual(actualVC.get(i), expectedVC.get(i), Collections.<String>emptyList(), Collections.<String>emptyList());
        }
    }

    @Test (dataProvider = "VCFdata")
    public void testUsingGenomicsDB(File[] inputs, File expected, File expectedDb, List<SimpleInterval> intervals, List<String> additionalArguments, String reference) throws IOException {
        Utils.resetRandomGenerator();
        //we need this for reproducibility because of random number generation in fixTooHighQD

        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(Arrays.asList(inputs), IntervalUtils.getSpanningInterval(intervals));
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        final File output = createTempFile("GnarlyGenotyper", ".vcf");
        final File outputDatabase = createTempFile("GnarlyGenotyper.annotationDatabase", ".vcf");

        runTool(genomicsDBUri, intervals, reference, output, outputDatabase, additionalArguments);

        final List<VariantContext> expectedVC = getVariantContexts(expected);
        final List<VariantContext> actualVC = getVariantContexts(output);
        assertEqualVariants(actualVC, expectedVC);
        if (expectedDb != null) {
            final List<VariantContext> expectedDB = getVariantContexts(expectedDb);
            final List<VariantContext> actualDB = getVariantContexts(outputDatabase);
            assertEqualVariants(actualDB, expectedDB);
        }
    }

    protected void runTool(String input, List<SimpleInterval> intervals, String reference, File output, File outputDatabase, List<String> additionalArguments) {

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
        .add("V", input)
        .add("output-db", outputDatabase.getAbsolutePath())
        .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");
        args.addOutput(output);
        intervals.forEach(args::addInterval);

        additionalArguments.forEach(args::addRaw);

        runCommandLine(args);
    }

    @Test
    public void testOnHailOutput() {
        final String input = getToolTestDataDir() + "hailOutput.chr20snippet.sites_only.vcf";
        final File output = createTempFile("GnarlyGenotyper", ".vcf");
        final File expected = new File(getToolTestDataDir() + "expected.hailOutput.chr20snippet.sites_only.vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .add("V", input)
                .add("L", "chr20:10000000-10030000")
                .add("only-output-calls-starting-in-intervals", true)
                .add("keep-all-sites", true)
                .addOutput(output)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");
        runCommandLine(args);

        //should have same variants as input with low QUAL variants marked with a monomorphic filter
        //QUAL column filled in, retains SB_TABLE, MQ_DP (deprecated), VarDP
        //does not expect ExcessHet because legacy input has neither genotypes nor GT counts
        try (final FeatureDataSource<VariantContext> actualVcs = new FeatureDataSource<>(output);
             final FeatureDataSource<VariantContext> expectedVcs = new FeatureDataSource<>(expected)) {
            GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                    (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                            Collections.emptyList(), Collections.emptyList()));
        }
    }
}