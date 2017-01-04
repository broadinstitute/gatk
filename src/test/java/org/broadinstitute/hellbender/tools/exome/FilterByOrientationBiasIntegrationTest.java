package org.broadinstitute.hellbender.tools.exome;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.OrientationBiasFilterConstants;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.OrientationBiasUtils;
import org.broadinstitute.hellbender.tools.picard.analysis.artifacts.Transition;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class FilterByOrientationBiasIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_RESOURCE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/orientationbiasvariantfilter/");
    public static final String smallM2Vcf = TEST_RESOURCE_DIR.getAbsolutePath() + "/small_m2.vcf";
    public static final String smallM2VcfMore = TEST_RESOURCE_DIR.getAbsolutePath() + "/small_m2_more_variants.vcf";
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
        final List<VariantContext> variantContexts = new ArrayList<>();
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(outputFile);
        for (final VariantContext vc : featureDataSource) {
            variantContexts.add(vc);
        }

        Assert.assertEquals(variantContexts.size(), 11);

        // Make sure that every entry has a orientation_bias filter in the genotype on the TUMOR sample if G/T or C/A
        for (final VariantContext vc: variantContexts) {
            final Genotype tumorGenotype = vc.getGenotype("TUMOR");
            Assert.assertTrue((tumorGenotype.getFilters() == null) || (tumorGenotype.getFilters().contains(OrientationBiasFilterConstants.IS_ORIENTATION_BIAS_CUT)) ||
                    !OrientationBiasUtils.isGenotypeInTransitionWithComplement(tumorGenotype, Transition.transitionOf('G', 'T')));

            final Genotype normalGenotype = vc.getGenotype("NORMAL");
            Assert.assertTrue((normalGenotype.getFilters() == null)
                    || normalGenotype.getFilters().equals(VCFConstants.UNFILTERED)
                    || normalGenotype.getFilters().equals(VCFConstants.PASSES_FILTERS_v4));
        }
    }
}
