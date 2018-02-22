package org.broadinstitute.hellbender.tools;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class AnnotatePairOrientationIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_BAM_DIR = new File("src/test/resources/large/mutect/dream_synthetic_bams/");

    final static String TEST_VCF = toolsTestDir + "/test_no_pair_orientation_info.vcf";
    final static String TEST_VCF_INDELS = toolsTestDir + "/test_no_pair_orientation_info_indels.vcf";
    final static String TEST_BAM_TUMOR = TEST_BAM_DIR.getAbsolutePath() + "/tumor_1.bam";
    final static String TEST_BAM_NORMAL = TEST_BAM_DIR.getAbsolutePath() + "/normal_1.bam";
    final static String TEST_BAM_TUMOR_INDELS = TEST_BAM_DIR.getAbsolutePath() + "/tumor_3.bam";
    final static String TEST_BAM_NORMAL_INDELS = TEST_BAM_DIR.getAbsolutePath() + "/normal_3.bam";

    // TODO: Test with multiallelics
    // TODO: Test with multiallelics and symbolic at the same time
    // TODO: Test with symbolic
    // TODO: Test with information missing from the VCF and make sure appropriate exception is thrown.
    // TODO: Test with more cutoff variables
    // TODO: Once above five TODOs are done (at least), AnnotatePairOrientation can be taken out of Experimental status.


    @Test
    public void testBasicIndels() throws IOException {
        final File outputFile = File.createTempFile("ob_indel_annotate_", ".vcf");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(TEST_VCF_INDELS);
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(TEST_BAM_TUMOR_INDELS);
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(TEST_BAM_NORMAL_INDELS);


        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        // Ground truth from manual review in IGV
        final String[][] gtF1R2F2R1 = {{"14,0", "21,0", "11,3", "13,4"},{"34,0","27,0", "10,12","14,11"},
                {"14,0", "14,0", "18,1", "17,3"},{"24,0","15,0", "19,7","22,2"}};

        Assert.assertTrue(outputFile.exists());
        final List<VariantContext> variantContexts = getVariantContextsFromFile(outputFile);

        assertOrientationAnnotationValues(variantContexts, gtF1R2F2R1, "G15512.prenormal.sorted",
                "IS3.snv.indel.sv");
    }

    /**
     *  Only tests SNVs
     * @throws IOException
     */
    @Test
    public void testBasicRun() throws IOException{
        final File outputFile = File.createTempFile("ob_annotate_", ".vcf");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(TEST_VCF);
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(TEST_BAM_TUMOR);
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(TEST_BAM_NORMAL);

        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        final List<VariantContext> variantContexts = getVariantContextsFromFile(outputFile);

        // Ground truth from manual review in IGV
        final String[][] gtF1R2F2R1 = {{"22,0", "11,0", "9,9", "8,4"},{"11,0","15,0", "11,8","9,10"}};

        assertOrientationAnnotationValues(variantContexts, gtF1R2F2R1, "synthetic.challenge.set1.normal",
                "tumor sample");
    }

    private List<VariantContext> getVariantContextsFromFile(File vcfFile) {
        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(vcfFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }
        return variantContexts;
    }

    private void assertOrientationAnnotationValues(final List<VariantContext> variantContexts, final String[][] gtF1R2F2R1,
                                                   final String normalSampleName, final String tumorSampleName) {
        final List<String> annotations = new ArrayList<>();
        annotations.add(GATKVCFConstants.F1R2_KEY);
        annotations.add(GATKVCFConstants.F2R1_KEY);

        for (int i = 0; i < variantContexts.size(); i++) {
            final VariantContext vc = variantContexts.get(i);
            final Genotype normalGenotype = vc.getGenotype(normalSampleName);
            Assert.assertTrue(normalGenotype.hasExtendedAttribute(GATKVCFConstants.F1R2_KEY));
            Assert.assertTrue(normalGenotype.hasExtendedAttribute(GATKVCFConstants.F2R1_KEY));

            for (int j = 0; j < annotations.size(); j ++) {
                final String annotation = annotations.get(j);
                final String normalF1r2 = normalGenotype.getExtendedAttribute(annotation).toString();
                Assert.assertEquals(normalF1r2, gtF1R2F2R1[i][j]);
            }

            final Genotype tumorGenotype = vc.getGenotype(tumorSampleName);
            Assert.assertTrue(tumorGenotype.hasExtendedAttribute(GATKVCFConstants.F1R2_KEY));
            Assert.assertTrue(tumorGenotype.hasExtendedAttribute(GATKVCFConstants.F2R1_KEY));
            for (int j = 0; j < annotations.size(); j ++) {
                final String annotation = annotations.get(j);
                final String tumorF1r2 = tumorGenotype.getExtendedAttribute(annotation).toString();
                Assert.assertNotNull(Utils.split(tumorF1r2, ","));
                Assert.assertEquals(Utils.split(tumorF1r2, ",").size(), 2);
                Assert.assertEquals(tumorF1r2, gtF1R2F2R1[i][j+annotations.size()]);
            }
        }
    }
}
