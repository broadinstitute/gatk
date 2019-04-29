package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBConstants;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
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
                // Simple Test, spanning deletions
               // {new File[]{getTestFile("sample1.vcf"), getTestFile("sample2.vcf"), getTestFile("sample3.vcf"), getTestFile("sample4.vcf"), getTestFile("sample5.vcf")},
               //         getTestFile("fiveSampleTest.vcf"), Arrays.asList(new SimpleInterval("chr20", 251370, 252000), new SimpleInterval("chr20", 253000, 255600)), Arrays.asList("--merge-input-intervals", "--only-output-calls-starting-in-intervals"), b38_reference_20_21},
               // {new File[]{getTestFile("sample1.vcf"), getTestFile("sample2.vcf"), getTestFile("sample3.vcf"), getTestFile("sample4.vcf"), getTestFile("sample5.vcf")},
               //         getTestFile("fiveSampleTest.vcf"), Arrays.asList(new SimpleInterval("chr20", 250865, 348163)), Arrays.asList("-stand-call-conf 10"), b38_reference_20_21},
                {new File[]{new File(getToolTestDataDir() + "/../variantutils/ReblockGVCF/expected.NA12878.AS.chr20snippet.reblocked.g.vcf"),
                        new File(getToolTestDataDir() + "/../variantutils/ReblockGVCF/expected.NA12892.AS.chr20snippet.reblocked.g.vcf")},
                        getTestFile("twoSampleAS.vcf"), Arrays.asList(new SimpleInterval("20")), NO_EXTRA_ARGS, b37_reference_20_21
                }
        };
    }


    @Test (dataProvider = "VCFdata")
    public void testUsingGenomicsDB(File[] inputs, File expected, List<SimpleInterval> intervals, List<String> additionalArguments, String reference) throws IOException {
        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(Arrays.asList(inputs), IntervalUtils.getSpanningInterval(intervals));
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        File output = runTool(genomicsDBUri, intervals, reference, additionalArguments);

        final List<VariantContext> expectedVC = getVariantContexts(expected);
        final List<VariantContext> actualVC = getVariantContexts(output);
        VariantContextTestUtils.assertEqualVariants(actualVC, expectedVC);
    }

    protected File runTool(String input, List<SimpleInterval> intervals, String reference, List<String> additionalArguments) {
        final File output = createTempFile("GnarlyGenotyper", ".vcf");
        final File outputDatabase = createTempFile("GnarlyGenotyper.annotationDatabase", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
        .addArgument("V", input)
        .addArgument("output-db", outputDatabase.getAbsolutePath());
        args.addOutput(output);
        intervals.forEach(args::addInterval);

        additionalArguments.forEach(args::add);

        runCommandLine(args);
        return output;
    }

    @Test
    public void testOnHailOutput() {
        final String input = getToolTestDataDir() + "hailTest.chr20snippet.sites_only.vcf";
        final File output = createTempFile("GnarlyGenotyper", ".vcf");
        final File expected = new File(getToolTestDataDir() + "expected.hailTest.chr20snippet.sites_only.vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg38Reference))
                .addArgument("V", input)
                .addArgument("L", "chr20:10000000-11000000")
                .addBooleanArgument("only-output-calls-starting-in-intervals", true)
                .addOutput(output);
        runCommandLine(args);

        try (final FeatureDataSource<VariantContext> actualVcs = new FeatureDataSource<>(output);
             final FeatureDataSource<VariantContext> expectedVcs = new FeatureDataSource<>(expected)) {
            GATKBaseTest.assertCondition(actualVcs, expectedVcs,
                    (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqual(a, e,
                            Collections.emptyList()));
        }
    }
}