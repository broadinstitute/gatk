package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.walkers.validation.AnnotateVcfWithBamDepth;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static org.testng.Assert.*;

/**
 * Created by tsato on 7/6/17.
 */
public class AnnotateVcfWithGCContentIntegrationTest extends CommandLineProgramTest {

    private static final String DREAM_BAMS_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final String DREAM_VCFS_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/dream/vcfs/";

    // test on the DREAM bam 1 and accompanying variants
    // depths verified manually in IGV
    @Test
    public void test() {
        final File vcf = new File(DREAM_VCFS_DIR, "sample_1.vcf");
        final File outputVcf = createTempFile("annotated", ".vcf");

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, b37_reference_20_21,
                "-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME, vcf.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputVcf.getAbsolutePath()
        };

        runCommandLine(arguments);

        final List<VariantContext> input = StreamSupport.stream(new FeatureDataSource<VariantContext>(vcf).spliterator(), false)
                .collect(Collectors.toList());
        final List<VariantContext> output = StreamSupport.stream(new FeatureDataSource<VariantContext>(outputVcf).spliterator(), false)
                .collect(Collectors.toList());

        Assert.assertEquals(input.size(), output.size());
    }

}