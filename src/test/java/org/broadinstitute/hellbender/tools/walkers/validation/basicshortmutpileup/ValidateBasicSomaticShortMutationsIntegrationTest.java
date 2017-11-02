package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.ValidateBasicSomaticShortMutations.*;

public class ValidateBasicSomaticShortMutationsIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DREAM_BAM_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final String TEST_DREAM_NORMAL_BAM = TEST_DREAM_BAM_DIR + "normal.bam";
    private static final String TEST_DREAM_TUMOR_BAM = TEST_DREAM_BAM_DIR + "tumor.bam";
    private static final String TEST_DREAM_VCF = toolsTestDir + "/walkers/validation/basicshortmutpileup/synthetic.challenge.set1.tumor-vs-synthetic.challenge.set1.normal-filtered.vcf";
    private static final String TEST_DREAM_NORMAL_BAM_INDEL_TEST = TEST_DREAM_BAM_DIR + "normal_3.bam";
    private static final String TEST_DREAM_TUMOR_BAM_INDEL_TEST = TEST_DREAM_BAM_DIR + "tumor_3.bam";
    private static final String TEST_DREAM_VCF_INDEL_TEST = toolsTestDir + "/walkers/validation/basicshortmutpileup/IS3.snv.indel.sv-vs-G15512.prenormal.sorted.vcf";
    private static final String REFERENCE = largeFileTestDir + "human_g1k_v37.20.21.fasta";

    @Test
    public void testBasic() {
        // This test is simply running the full tool and making sure that there are no serious errors.
        //  No variants should validate, since the validation bam is not the same one used for calling.
        final File outputFile = IOUtils.createTempFile("basicTest", ".txt");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ValidateBasicSomaticShortMutations.SAMPLE_NAME_DISCOVERY_VCF);
        arguments.add("synthetic.challenge.set1.tumor");
        arguments.add("-" + ValidateBasicSomaticShortMutations.SAMPLE_NAME_VALIDATION_CASE);
        arguments.add("synthetic.challenge.set1.tumor");
        arguments.add("-" + ValidateBasicSomaticShortMutations.SAMPLE_NAME_VALIDATION_CONTROL);
        arguments.add("synthetic.challenge.set1.normal");

        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(TEST_DREAM_VCF);

        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(TEST_DREAM_TUMOR_BAM);
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(TEST_DREAM_NORMAL_BAM);

        arguments.add("-L");
        arguments.add(TEST_DREAM_VCF);

        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REFERENCE);

        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        arguments.add("--verbosity");
        arguments.add("INFO");
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<SimpleAnnotatedGenomicRegion> variantValidationResults =
                SimpleAnnotatedGenomicRegion.readAnnotatedRegions(outputFile, new HashSet<>(Arrays.asList(ValidateBasicSomaticShortMutations.headers)));

        Assert.assertEquals(variantValidationResults.size(), 2);

        Assert.assertEquals(variantValidationResults.get(0).getInterval(), new SimpleInterval("20", 10022820, 10022820));
        // The variant in the VCF is A>C, the bam file has A>T, so alt count should be zero.
        Assert.assertEquals(variantValidationResults.get(0).getAnnotations().get(VALIDATION_ALT_COVERAGE), "0");
        Assert.assertEquals(variantValidationResults.get(0).getAnnotations().get(VALIDATION_REF_COVERAGE), "18");
        Assert.assertEquals(variantValidationResults.get(0).getAnnotations().get(IS_NOT_NOISE), "false");
        Assert.assertEquals(variantValidationResults.get(0).getAnnotations().get(DISCOVERY_VCF_FILTER), "");

        Assert.assertEquals(variantValidationResults.get(1).getInterval(), new SimpleInterval("20", 10080550, 10080550));
        Assert.assertEquals(variantValidationResults.get(1).getAnnotations().get(VALIDATION_ALT_COVERAGE), "0");
        Assert.assertEquals(variantValidationResults.get(1).getAnnotations().get(VALIDATION_REF_COVERAGE), "16");
        Assert.assertEquals(variantValidationResults.get(1).getAnnotations().get(IS_NOT_NOISE), "false");
        Assert.assertEquals(variantValidationResults.get(1).getAnnotations().get(DISCOVERY_VCF_FILTER), "artifact_in_normal;germline_risk;t_lod");
    }

    @Test
    public void testIndelCase() {
        // This test is simply running the full tool and making sure that there are no serious errors.
        //  All variants should validate (or be unpowered)

        final File outputFile = IOUtils.createTempFile("basicTest", ".txt");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ValidateBasicSomaticShortMutations.SAMPLE_NAME_DISCOVERY_VCF);
        arguments.add("IS3.snv.indel.sv");
        arguments.add("-" + ValidateBasicSomaticShortMutations.SAMPLE_NAME_VALIDATION_CASE);
        arguments.add("IS3.snv.indel.sv");
        arguments.add("-" + ValidateBasicSomaticShortMutations.SAMPLE_NAME_VALIDATION_CONTROL);
        arguments.add("G15512.prenormal.sorted");

        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(TEST_DREAM_VCF_INDEL_TEST);

        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(TEST_DREAM_NORMAL_BAM_INDEL_TEST);
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(TEST_DREAM_TUMOR_BAM_INDEL_TEST);

        arguments.add("-L");
        arguments.add(TEST_DREAM_VCF_INDEL_TEST);

        arguments.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
        arguments.add(REFERENCE);
        arguments.add("-" + ValidateBasicSomaticShortMutations.CUTOFF_LONG_NAME);
        arguments.add(String.valueOf(ValidateBasicSomaticShortMutations.DEFAULT_MIN_BQ_CUTOFF));

        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        arguments.add("--verbosity");
        arguments.add("INFO");
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<SimpleAnnotatedGenomicRegion> variantValidationResults =
                SimpleAnnotatedGenomicRegion.readAnnotatedRegions(outputFile, new HashSet<>(Arrays.asList(ValidateBasicSomaticShortMutations.headers)));

        Assert.assertEquals(variantValidationResults.size(), 336);


        // DEL: 20	1330646	.	CCTTGGCTTATTCCA	C
        // INS1: 20	3076247	.	AT	ATT
        // INS2: 20	3076299	.	A	AAGAAGCATGC

        assertValidationResult(variantValidationResults, new SimpleInterval("20", 1330646, 1330646),
                "CCTTGGCTTATTCCA", "C", 7, 21, 10,
                26);
        assertValidationResult(variantValidationResults, new SimpleInterval("20", 3076247, 3076247),
                "AT", "ATT", 4, 33, 4,
                25);
        assertValidationResult(variantValidationResults, new SimpleInterval("20", 3076299, 3076299),
                "A", "AAGAAGCATGC", 9, 41, 12,
                37);  // One read has low BQ in the insertion, so 9, instead of 10
    }

    private void assertValidationResult(final List<SimpleAnnotatedGenomicRegion> variantValidationResults,
                                        final SimpleInterval firstBaseInVariant, final String refString,
                                        final String altString, int gtValidationAltCoverage, int gtValidationRefCoverage,
                                        int gtDiscoveryAltCoverage, int gtDiscoveryRefCoverage) {
        final OverlapDetector<SimpleAnnotatedGenomicRegion> overlapDetector = OverlapDetector.create(variantValidationResults);
        final SimpleAnnotatedGenomicRegion indel = overlapDetector.getOverlaps(firstBaseInVariant).iterator().next();
        final SortedMap<String, String> indelAnnotations = indel.getAnnotations();

        Assert.assertEquals(indelAnnotations.get(ValidateBasicSomaticShortMutations.REF), refString);
        Assert.assertEquals(indelAnnotations.get(ValidateBasicSomaticShortMutations.ALT), altString);

        final int validationAltCoverage = Integer.parseInt(indelAnnotations.get(ValidateBasicSomaticShortMutations.VALIDATION_ALT_COVERAGE));
        final int validationRefCoverage = Integer.parseInt(indelAnnotations.get(ValidateBasicSomaticShortMutations.VALIDATION_REF_COVERAGE));
        final int discoveryAltCoverage = Integer.parseInt(indelAnnotations.get(ValidateBasicSomaticShortMutations.DISCOVERY_ALT_COVERAGE));
        final int discoveryRefCoverage = Integer.parseInt(indelAnnotations.get(ValidateBasicSomaticShortMutations.DISCOVERY_REF_COVERAGE));
        final double power = Double.parseDouble(indelAnnotations.get(ValidateBasicSomaticShortMutations.POWER));
        final int minCount = Integer.parseInt(indelAnnotations.get(ValidateBasicSomaticShortMutations.MIN_VAL_COUNT));

        Assert.assertEquals(validationAltCoverage, gtValidationAltCoverage);
        Assert.assertEquals(validationRefCoverage, gtValidationRefCoverage);
        Assert.assertEquals(discoveryAltCoverage, gtDiscoveryAltCoverage);
        Assert.assertEquals(discoveryRefCoverage, gtDiscoveryRefCoverage);
        Assert.assertEquals(power,
                PowerCalculationUtils.calculatePower(validationAltCoverage + validationRefCoverage,
                        discoveryAltCoverage, discoveryAltCoverage + discoveryRefCoverage,
                        minCount));
    }
}
