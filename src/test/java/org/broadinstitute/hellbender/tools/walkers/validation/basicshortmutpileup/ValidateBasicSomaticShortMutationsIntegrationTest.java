package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

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
        final File outputFile = IOUtils.createTempFile("basicTest", ".seg");
        final File summaryFile = IOUtils.createTempFile("summary", ".txt");
        final File annotatedVcf = IOUtils.createTempFile("annotated", ".vcf");
        final List<String> arguments = new ArrayList<>();
        arguments.add("--" + ValidateBasicSomaticShortMutations.SAMPLE_NAME_DISCOVERY_VCF_LONG_NAME);
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
        arguments.add("-" + Concordance.SUMMARY_LONG_NAME);
        arguments.add(summaryFile.getAbsolutePath());
        arguments.add("--" + ValidateBasicSomaticShortMutations.ANNOTATED_VCF_LONG_NAME);
        arguments.add(annotatedVcf.getAbsolutePath());
        arguments.add("--verbosity");
        arguments.add("INFO");
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(summaryFile.exists());

        final List<BasicValidationResult> variantValidationResults = BasicValidationResult.read(new GATKPath(outputFile.getAbsolutePath()));

        //final List<AnnotatedInterval> variantValidationResults =
        //        AnnotatedIntervalCollection.create(outputFile.toPath(), new HashSet<>(Arrays.asList(ValidateBasicSomaticShortMutations.headers))).getRecords();

        Assert.assertEquals(variantValidationResults.size(), 2);

        Assert.assertEquals(variantValidationResults.get(0).getInterval(), new SimpleInterval("20", 10022820, 10022820));
        // The variant in the VCF is A>C, the bam file has A>T, so alt count should be zero.
        Assert.assertEquals(variantValidationResults.get(0).getValidationAltCount(), 0);
        Assert.assertEquals(variantValidationResults.get(0).getValidationRefCount(), 18);
        Assert.assertFalse(variantValidationResults.get(0).isOutOfNoiseFloor());
        Assert.assertEquals(variantValidationResults.get(0).getFilters(), "");
        Assert.assertEquals(variantValidationResults.get(0).getNumAltSupportingReadsInNormal(), 0);

        Assert.assertEquals(variantValidationResults.get(1).getInterval(), new SimpleInterval("20", 10080550, 10080550));
        Assert.assertEquals(variantValidationResults.get(1).getValidationAltCount(), 0);
        Assert.assertEquals(variantValidationResults.get(1).getValidationRefCount(), 16);
        Assert.assertFalse(variantValidationResults.get(1).isOutOfNoiseFloor());
        Assert.assertEquals(variantValidationResults.get(1).getFilters(), "artifact_in_normal;germline_risk;t_lod");
        Assert.assertEquals(variantValidationResults.get(1).getNumAltSupportingReadsInNormal(), 0);
    }

    @Test
    public void testIndelCase() {
        // This test is simply running the full tool and making sure that there are no serious errors.
        //  All variants should validate (or be unpowered)

        final File outputFile = IOUtils.createTempFile("basicTest", ".seg");
        final File summaryFile = IOUtils.createTempFile("summary", ".txt");

        final List<String> arguments = new ArrayList<>();
        arguments.add("--" + ValidateBasicSomaticShortMutations.SAMPLE_NAME_DISCOVERY_VCF_LONG_NAME);
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
        arguments.add("-" + Concordance.SUMMARY_LONG_NAME);
        arguments.add(summaryFile.getAbsolutePath());

        arguments.add("--verbosity");
        arguments.add("INFO");
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(summaryFile.exists());

        final List<BasicValidationResult> variantValidationResults = BasicValidationResult.read(new GATKPath(outputFile.getAbsolutePath()));

        Assert.assertEquals(variantValidationResults.size(), 336);

        // DEL: 20	1330646	.	CCTTGGCTTATTCCA	C
        // INS1: 20	3076247	.	AT	ATT
        // INS2: 20	3076299	.	A	AAGAAGCATGC

        assertValidationResult(variantValidationResults, new SimpleInterval("20", 1330646, 1330646),
                "CCTTGGCTTATTCCA", "C", 7, 21, 10,
                26, 0);

        // TODO: This next one actually has an incorrect gtNumAltReadsInValidationNormal, since the validator thinks of this as AT<insT>.  There are supporting reads in the validation normal for A<insT>T.  See https://github.com/broadinstitute/gatk/issues/5061
        assertValidationResult(variantValidationResults, new SimpleInterval("20", 3076247, 3076247),
                "AT", "ATT", 4, 33, 4,
                25,0);
        assertValidationResult(variantValidationResults, new SimpleInterval("20", 3076299, 3076299),
                "A", "AAGAAGCATGC", 9, 41, 12,
                37, 0);  // One read has low BQ in the insertion, so 9, instead of 10
    }

    private void assertValidationResult(final List<BasicValidationResult> variantValidationResults,
                                        final SimpleInterval firstBaseInVariant, final String refString,
                                        final String altString, int gtValidationAltCoverage, int gtValidationRefCoverage,
                                        int gtDiscoveryAltCoverage, int gtDiscoveryRefCoverage, int gtNumAltReadsInValidationNormal) {
        final OverlapDetector<BasicValidationResult> overlapDetector = OverlapDetector.create(variantValidationResults);
        final BasicValidationResult indel = overlapDetector.getOverlaps(firstBaseInVariant).iterator().next();

        Assert.assertEquals(indel.getReference().getBaseString(), refString);
        Assert.assertTrue(indel.getReference().isReference());
        Assert.assertTrue(indel.getAlternate().isNonReference());
        Assert.assertEquals(indel.getAlternate().getBaseString(), altString);

        final int validationAltCoverage = indel.getValidationAltCount();
        final int validationRefCoverage = indel.getValidationRefCount();
        final int discoveryAltCoverage = indel.getDiscoveryAltCount();
        final int discoveryRefCoverage = indel.getDiscoveryRefCount();
        final double power = indel.getPower();
        final int minCount = indel.getMinValidationReadCount();
        final long numSupportingAltReadsInValidationNormal = indel.getNumAltSupportingReadsInNormal();

        Assert.assertEquals(validationAltCoverage, gtValidationAltCoverage);
        Assert.assertEquals(validationRefCoverage, gtValidationRefCoverage);
        Assert.assertEquals(discoveryAltCoverage, gtDiscoveryAltCoverage);
        Assert.assertEquals(discoveryRefCoverage, gtDiscoveryRefCoverage);
        Assert.assertEquals(power,
                PowerCalculationUtils.calculatePower(validationAltCoverage + validationRefCoverage,
                        discoveryAltCoverage, discoveryAltCoverage + discoveryRefCoverage,
                        minCount));
        Assert.assertEquals(numSupportingAltReadsInValidationNormal, gtNumAltReadsInValidationNormal, "Discordance in alt count in validation normal at " + firstBaseInVariant);
    }
}
