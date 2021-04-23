package org.broadinstitute.hellbender.tools;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.SparkTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerIntegrationTest;

@Test(groups = {"variantcalling"})
// Selected tests copied from HaplotypeCallerIntegrationTest
public class HaplotypeCallerSparkIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=HaplotypeCallerIntegrationTest"
    // to update all of the exact-match tests at once. After you do this, you should look at the
    // diffs in the new expected outputs in git to confirm that they are consistent with expectations.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    public static final String TEST_FILES_DIR = toolsTestDir + "haplotypecaller/";

    /*
     * Make sure that someone didn't leave the UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS toggle turned on
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @DataProvider(name = "HaplotypeCallerTestInputs")
    public Object[][] getHaplotypCallerTestInputs() {
        return new Object[][]{
                {NA12878_20_21_WGS_bam, b37_reference_20_21},
                // Uncomment when bai indexed CRAM is supported in Disq
                //{NA12878_20_21_WGS_cram, b37_reference_20_21}
        };
    }

    /*
     * Test that in VCF mode we're consistent with past GATK4 results
     */
    @Test(dataProvider = "HaplotypeCallerTestInputs")
    public void testVCFModeIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testVCFMode.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "--strict" // to match walker version
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    @Test(dataProvider = "HaplotypeCallerTestInputs")
    public void testNonStrictVCFModeIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testVCFMode.gatk4.nonstrict.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    /*
    * Test that in VCF mode we're >= 99% concordant with GATK3.8 results
    */
    @Test
    public void testVCFModeIsConcordantWithGATK3_8Results() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConcordantWithGATK3Results", ".vcf");
        //Created by running:
        // java -jar gatk.3.8-4-g7b0250253f.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // --out expected.testVCFMode.gatk3.8-4-g7b0250253.vcf -G StandardHC -G Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk3.8-4-g7b0250253.vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        final double concordance = HaplotypeCallerIntegrationTest.calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in VCF mode is < 99% (" +  concordance + ")");
    }

    /**
     * Test that in VCF mode we're >= 99% concordant with GATK3.8 results
     * THIS TEST explodes with an exception because Allele-Specific annotations are not supported in vcf mode yet.
     * It's included to parallel the matching (also exploding) test for the non-spark HaplotypeCaller
     * {@link org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerIntegrationTest#testVCFModeIsConcordantWithGATK3_8ResultsAlleleSpecificAnnotations(String, String)}}
     */
    @Test(expectedExceptions = UserException.class)
    public void testVCFModeIsConcordantWithGATK3_8ResultsAlleleSpecificAnnotations() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConcordantWithGATK3.8ResultsAlleleSpecificAnnotations", ".vcf");

        //Created by running
        //java -jar gatk.3.8-4-g7b0250253f.jar  -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // --out expected.testVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.vcf -G StandardHC -G Standard -G AS_Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "AS_StandardAnnotation",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        final double concordance = HaplotypeCallerIntegrationTest.calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in AS VCF mode is < 99% (" +  concordance + ")");
    }

    /*
   * Test that in GVCF mode we're >= 99% concordant with GATK3 results
   */
    @Test(enabled=false) //disabled after reference confidence change in #5172
    public void testGVCFModeIsConcordantWithGATK3_8Results() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConcordantWithGATK3Results", ".g.vcf");
        //Created by running:
        //java -jar  gatk.3.8-4-g7b0250253f.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // -ERC GVCF --out expected.testGVCFMode.3.8-4-g7b0250253f.g.vcf -G StandardHC -G Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testGVCFMode.3.8-4-g7b0250253f.g.vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, ReferenceConfidenceMode.GVCF.toString(),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        final double concordance = HaplotypeCallerIntegrationTest.calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in GVCF mode is < 99% (" +  concordance + ")");
    }

    @DataProvider
    public static Object[][] brokenGVCFCases() {
        return new Object[][]{
                {".g.bcf"},
                {".g.bcf.gz"}
        };
    }

    @Test(dataProvider = "brokenGVCFCases", expectedExceptions = UserException.UnimplementedFeature.class)
    public void testBrokenGVCFConfigurationsAreDisallowed(String extension) {
        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-O", createTempFile("testGVCF_GZ_throw_exception", extension).getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, ReferenceConfidenceMode.GVCF.toString(),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
        };

        runCommandLine(args);
    }

    @DataProvider
    public static Object[][] gvcfCases() {
        return new Object[][]{
                {".g.vcf"},
                {".g.vcf.gz"}
        };
    }

    @Test(dataProvider = "gvcfCases", enabled=false) //disabled after reference confidence change in #5172
    public void testGVCFModeIsConcordantWithGATK3_8AlelleSpecificResults(String extension) throws Exception {
        Utils.resetRandomGenerator();
        final File output = createTempFile("testGVCFModeIsConcordantWithGATK3_8AlelleSpecificResults", extension);

        //Created by running:
        // java -jar gatk.3.8-4-g7b0250253f.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // -ERC GVCF --out expected.testGVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.g.vcf -G StandardHC -G Standard -G AS_Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testGVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.g.vcf");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "AS_StandardAnnotation",
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, ReferenceConfidenceMode.GVCF.toString(),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        final double concordance = HaplotypeCallerIntegrationTest.calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in AS GVCF mode is < 99% (" +  concordance + ")");
    }

    @Test
    public void testGenotypeCalculationArgumentCollectionIsSerializable() {
        final GenotypeCalculationArgumentCollection args = new GenotypeCalculationArgumentCollection();
        SparkTestUtils.roundTripInKryo(args, GenotypeCalculationArgumentCollection.class, SparkContextFactory.getTestSparkContext().getConf());

    }

    @Test
    public void testHaplotypeCallerArgsAreSerializable() {
        final HaplotypeCallerArgumentCollection args = new HaplotypeCallerArgumentCollection();
        SparkTestUtils.roundTripInKryo(args, HaplotypeCallerArgumentCollection.class, SparkContextFactory.getTestSparkContext().getConf());
    }


    @Test
    public void testReferenceMultiSourceIsSerializable() {
        final ReferenceMultiSparkSource args = new ReferenceMultiSparkSource(new GATKPath(GATKBaseTest.b37_reference_20_21), ReferenceWindowFunctions.IDENTITY_FUNCTION);
        SparkTestUtils.roundTripInKryo(args, ReferenceMultiSparkSource.class, SparkContextFactory.getTestSparkContext().getConf());
    }


    @Test
    public void testBroadcastHcArgs() {
        Broadcast<HaplotypeCallerArgumentCollection> broadcast = SparkContextFactory.getTestSparkContext().broadcast(new HaplotypeCallerArgumentCollection());
        broadcast.getValue();
    }

    @Test
    public void testFastGenotypeIsSerializable() {
        Genotype genotype = GenotypeBuilder.create("sample1", Collections.nCopies(2, Allele.create("C", false)));
        SparkTestUtils.roundTripInKryo(genotype, genotype.getClass(), SparkContextFactory.getTestSparkContext().getConf());
    }

    @Test
    public void testAllelesAreSerializable() {
        Allele a = Allele.create("A");
        SparkTestUtils.roundTripInKryo(a, a.getClass(), SparkContextFactory.getTestSparkContext().getConf());
        SparkTestUtils.roundTripInKryo(Allele.NO_CALL, Allele.class, SparkContextFactory.getTestSparkContext().getConf());
    }
}
