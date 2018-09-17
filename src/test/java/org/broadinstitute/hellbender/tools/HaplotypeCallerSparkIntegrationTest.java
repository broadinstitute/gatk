package org.broadinstitute.hellbender.tools;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.SparkTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerIntegrationTest;

@Test(groups = {"variantcalling"})
public class HaplotypeCallerSparkIntegrationTest extends CommandLineProgramTest {

    private  static final String TEST_FILES_DIR = toolsTestDir + "haplotypecaller/";

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
                "-R", b37_2bit_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
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
     * {@link org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerIntegrationTest#testVCFModeIsConcordantWithGATK3_8ResultsAlleleSpecificAnnotations()}
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
                "-R", b37_2bit_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "AS_StandardAnnotation",
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
                "-R", b37_2bit_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-ERC", "GVCF",
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
                "-R", b37_2bit_reference_20_21,
                "-O", createTempFile("testGVCF_GZ_throw_exception", extension).getAbsolutePath(),
                "-ERC", "GVCF",
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
                "-R", b37_2bit_reference_20_21,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "AS_StandardAnnotation",
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING"
        };

        runCommandLine(args);

        final double concordance = HaplotypeCallerIntegrationTest.calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in AS GVCF mode is < 99% (" +  concordance + ")");
    }

    @Test
    public void testReferenceAdapterIsSerializable() throws IOException {
        final ReferenceMultiSparkSource referenceMultiSource = new ReferenceMultiSparkSource(b37_2bit_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
        SparkTestUtils.roundTripInKryo(referenceMultiSource, ReferenceMultiSparkSource.class, SparkContextFactory.getTestSparkContext().getConf());
        final HaplotypeCallerSpark.ReferenceMultiSourceAdapter adapter = new HaplotypeCallerSpark.ReferenceMultiSourceAdapter(referenceMultiSource);
        SparkTestUtils.roundTripInKryo(adapter, HaplotypeCallerSpark.ReferenceMultiSourceAdapter.class, SparkContextFactory.getTestSparkContext().getConf());

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
        final ReferenceMultiSparkSource args = new ReferenceMultiSparkSource(GATKBaseTest.b37_2bit_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
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
