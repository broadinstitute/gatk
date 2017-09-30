package org.broadinstitute.hellbender.tools.exome;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.OrientationBiasFilterConstants;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.OrientationBiasUtils;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.OrientationSampleTransitionSummary;
import org.broadinstitute.hellbender.utils.artifacts.Transition;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class FilterByOrientationBiasIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_RESOURCE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/orientationbiasvariantfilter/");
    public static final String emptyVcf = TEST_RESOURCE_DIR.getAbsolutePath() + "/empty.vcf";
    public static final String emptyVcfNoSamples = TEST_RESOURCE_DIR.getAbsolutePath() + "/empty_and_no_samples.vcf";
    public static final String smallM2VcfMore = TEST_RESOURCE_DIR.getAbsolutePath() + "/small_m2_more_variants.vcf";
    public static final String smallHighDiploid = TEST_RESOURCE_DIR.getAbsolutePath() + "/high_ploidy.vcf";
    public static final String smallMA = TEST_RESOURCE_DIR.getAbsolutePath() + "/m2_multiallelic.vcf";
    public static final String nullADField = TEST_RESOURCE_DIR.getAbsolutePath() + "/null_AD_field.vcf";
    public static final String preAdapterQFile = TEST_RESOURCE_DIR.getAbsolutePath() + "/SAMPLE9.pre_adapter_detail_metrics";

    @Test
    public void testRun() throws IOException {
        final File outputFile = File.createTempFile("ob_", ".vcf");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + FilterByOrientationBias.PRE_ADAPTER_METRICS_DETAIL_FILE_SHORT_NAME);
        arguments.add(preAdapterQFile);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(smallM2VcfMore);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        final File summaryFile = new File(outputFile.getAbsolutePath() + FilterByOrientationBias.SUMMARY_FILE_SUFFIX);
        Assert.assertTrue(summaryFile.exists());
        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        Assert.assertEquals(variantContexts.size(), 11);
        Assert.assertTrue(FileUtils.sizeOf(outputFile) > 0);
        Assert.assertTrue(FileUtils.sizeOf(summaryFile) > 0);

        boolean is_variant_context_tested = false;

        // Make sure that every entry has a orientation_bias filter in the genotype on the TUMOR sample if G/T or C/A.
        //  Also, make sure that the variant context has the filter as well.  Not just the genotypes.
        for (final VariantContext vc: variantContexts) {
            final Genotype tumorGenotype = vc.getGenotype("TUMOR");
            Assert.assertTrue((tumorGenotype.getFilters() == null) || (tumorGenotype.getFilters().contains(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT)) ||
                    !OrientationBiasUtils.isGenotypeInTransitionWithComplement(tumorGenotype, Transition.transitionOf('G', 'T')));

            // If we see a filtered genotype, make sure the variant context was filtered as well.
            if ((tumorGenotype.getFilters() != null) && (tumorGenotype.getFilters().contains(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT))) {
                Assert.assertTrue(vc.getFilters().contains(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT));
                is_variant_context_tested = true;
            }

            final Genotype normalGenotype = vc.getGenotype("NORMAL");
            Assert.assertTrue((normalGenotype.getFilters() == null)
                    || normalGenotype.getFilters().equals(VCFConstants.UNFILTERED)
                    || normalGenotype.getFilters().equals(VCFConstants.PASSES_FILTERS_v4));
        }

        Assert.assertTrue(is_variant_context_tested, "Unit test may be broken.  Should have tested that variant context contained filter as well as genotype fields.");

        final List<OrientationSampleTransitionSummary> summaries = OrientationBiasUtils.readOrientationBiasSummaryTable(summaryFile);
        Assert.assertEquals(summaries.size(), 2);
        Assert.assertEquals(summaries.stream().filter(s -> s.getSample().equals("NORMAL")).count(), 1);
        Assert.assertEquals(summaries.stream().filter(s -> s.getSample().equals("TUMOR")).count(), 1);
        Assert.assertEquals(summaries.stream().filter(s -> s.getSample().equals("NORMAL")).map(s -> s.getArtifactMode()).filter(am -> am.equals(Transition.GtoT)).count(), 1);
        Assert.assertEquals(summaries.stream().filter(s -> s.getSample().equals("NORMAL")).map(s -> s.getArtifactModeComplement()).filter(am -> am.equals(Transition.CtoA)).count(), 1);
        Assert.assertEquals(summaries.stream().filter(s -> s.getSample().equals("TUMOR")).map(s -> s.getArtifactMode()).filter(am -> am.equals(Transition.GtoT)).count(), 1);
        Assert.assertEquals(summaries.stream().filter(s -> s.getSample().equals("TUMOR")).map(s -> s.getArtifactModeComplement()).filter(am -> am.equals(Transition.CtoA)).count(), 1);
        Assert.assertEquals(summaries.stream().filter(s -> s.getArtifactModeComplement().equals(s.getArtifactMode().complement())).count(), summaries.size());
        Assert.assertEquals(summaries.stream().filter(s -> s.getSample().equals("TUMOR")).map(s -> s.getArtifactModeComplement()).filter(am -> am.equals(Transition.CtoA)).count(), 1);
        Assert.assertEquals(summaries.stream().mapToLong(s -> s.getNumArtifactModeFiltered()).sum(), 4);

    }

    @Test
    public void testNullADField() throws IOException {
        final File outputFile = File.createTempFile("ob_", ".vcf");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + FilterByOrientationBias.PRE_ADAPTER_METRICS_DETAIL_FILE_SHORT_NAME);
        arguments.add(preAdapterQFile);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(nullADField);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        //Make sure we don't get an NPE
        Assert.assertTrue(outputFile.exists());
    }

    @Test
    public void testNoVariantsRun() throws IOException {
        final File outputFile = File.createTempFile("ob_no_variants", ".vcf");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + FilterByOrientationBias.PRE_ADAPTER_METRICS_DETAIL_FILE_SHORT_NAME);
        arguments.add(preAdapterQFile);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(emptyVcf);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        Assert.assertEquals(variantContexts.size(), 0);
        final File summaryFile = new File(outputFile.getAbsolutePath() + FilterByOrientationBias.SUMMARY_FILE_SUFFIX);
        Assert.assertTrue(summaryFile.exists());
        Assert.assertTrue(FileUtils.sizeOf(outputFile) > 0);
        Assert.assertTrue(FileUtils.sizeOf(summaryFile) > 0);

        final List<OrientationSampleTransitionSummary> summaries = OrientationBiasUtils.readOrientationBiasSummaryTable(summaryFile);
        Assert.assertEquals(summaries.size(), 0);
    }

    @Test
    public void testNoVariantsNoSamplesRun() throws IOException {
        final File outputFile = File.createTempFile("ob_no_variants_no_sample", ".vcf");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + FilterByOrientationBias.PRE_ADAPTER_METRICS_DETAIL_FILE_SHORT_NAME);
        arguments.add(preAdapterQFile);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(emptyVcfNoSamples);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        Assert.assertEquals(variantContexts.size(), 0);
        Assert.assertTrue(FileUtils.sizeOf(outputFile) > 0);
        final File summaryFile = new File(outputFile.getAbsolutePath() + FilterByOrientationBias.SUMMARY_FILE_SUFFIX);
        Assert.assertTrue(summaryFile.exists());
        Assert.assertTrue(FileUtils.sizeOf(summaryFile) > 0);

        final List<OrientationSampleTransitionSummary> summaries = OrientationBiasUtils.readOrientationBiasSummaryTable(summaryFile);
        Assert.assertEquals(summaries.size(), 0);
    }

    @Test
    public void testHighPloidyRun() throws IOException {
        final File outputFile = File.createTempFile("ob_high_ploidy", ".vcf");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + FilterByOrientationBias.PRE_ADAPTER_METRICS_DETAIL_FILE_SHORT_NAME);
        arguments.add(preAdapterQFile);
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(smallHighDiploid);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        Assert.assertEquals(variantContexts.size(), 1);
    }

    @Test
    public void testMultiallelic() throws IOException {
        final File outputFile = File.createTempFile("ob_ma", ".vcf");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + FilterByOrientationBias.PRE_ADAPTER_METRICS_DETAIL_FILE_SHORT_NAME);
        arguments.add(preAdapterQFile);
        arguments.add("-" + FilterByOrientationBias.ARTIFACT_MODES_SHORT_NAME);
        arguments.add("C/T");
        arguments.add("-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME);
        arguments.add(smallMA);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());
        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        // It is important to remember that the filter only looks at the first alternate allele.
        Assert.assertEquals(variantContexts.size(), 5);

        // The first variant should have a null in the been OB filter annotation.
        Assert.assertTrue(variantContexts.get(0).getGenotype(0).getExtendedAttribute("OBAMRC").toString().equals("false"));
        Assert.assertTrue(variantContexts.get(0).getGenotype(0).getExtendedAttribute("OBQRC").toString().equals("100.00"));

        // The second variant should not be filtered
        Assert.assertTrue(variantContexts.get(1).getGenotype(0).getExtendedAttribute("OBAMRC").toString().equals("true"));
        Assert.assertTrue(variantContexts.get(2).getGenotype(0).getExtendedAttribute("GT") == null);
        Assert.assertTrue(Double.parseDouble(variantContexts.get(1).getGenotype(0).getExtendedAttribute("OBQRC").toString()) < 100.0);

        // The third variant should be filtered
        Assert.assertTrue(variantContexts.get(2).getGenotype(0).getExtendedAttribute("OBAMRC").toString().equals("true"));
        Assert.assertTrue(Double.parseDouble(variantContexts.get(2).getGenotype(0).getExtendedAttribute("OBQRC").toString()) < 100.0);

        // fourth is in artifact mode
        Assert.assertTrue(variantContexts.get(3).getGenotype(0).getExtendedAttribute("OBAMRC").toString().equals("false"));
        Assert.assertTrue(variantContexts.get(3).getGenotype(0).getExtendedAttribute("OBAM").toString().equals("true"));
        Assert.assertTrue(Double.parseDouble(variantContexts.get(3).getGenotype(0).getExtendedAttribute("OBQ").toString()) < 100.0);

        Assert.assertTrue(variantContexts.get(4).getGenotype(0).getExtendedAttribute("OBAM").toString().equals("false"));
        Assert.assertTrue(variantContexts.get(4).getGenotype(0).getExtendedAttribute("OBQ").toString().equals("100.00"));

    }
}
