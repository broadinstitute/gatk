package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Created by davidben on 1/31/17.
 */
public class AnnotateVcfWithBamDepthIntegrationTest extends CommandLineProgramTest {
    private static final String DREAM_BAMS_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final String DREAM_VCFS_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/dream/vcfs/";

    // test on the DREAM bam 1 and accompanying variants
    // depths verified manually in IGV
    @Test
    public void test() {
        final File bam = new File(DREAM_BAMS_DIR, "tumor_1.bam");
        final File vcf = new File(DREAM_VCFS_DIR, "sample_1.vcf");
        final File outputVcf = createTempFile("annotated", ".vcf");

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME, vcf.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, bam.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputVcf.getAbsolutePath()
        };

        runCommandLine(arguments);

        final List<VariantContext> input = StreamSupport.stream(new FeatureDataSource<VariantContext>(vcf).spliterator(), false)
                .collect(Collectors.toList());
        final List<VariantContext> output = StreamSupport.stream(new FeatureDataSource<VariantContext>(outputVcf).spliterator(), false)
                .collect(Collectors.toList());

        Assert.assertEquals(input.size(), output.size());

        final List<String> inputKeys = input.stream().map(vc -> keyForVariant(vc)).collect(Collectors.toList());
        final List<String> outputKeys = output.stream().map(vc -> keyForVariant(vc)).collect(Collectors.toList());

        Assert.assertEquals(inputKeys, outputKeys);

        final List<Integer> bamDepths = output.stream()
                .map(vc -> vc.getAttributeAsInt(AnnotateVcfWithBamDepth.POOLED_BAM_DEPTH_ANNOTATION_NAME, -1))
                .collect(Collectors.toList());

        final List<Integer> firstSeveralDepthsFromIGV = Arrays.asList(33, 39, 19, 35, 25, 27);

        Assert.assertEquals(bamDepths.subList(0, firstSeveralDepthsFromIGV.size()), firstSeveralDepthsFromIGV);
    }

    private static String keyForVariant( final VariantContext variant ) {
        return String.format("%s:%d-%d %s", variant.getContig(), variant.getStart(), variant.getEnd(), variant.getAlleles());
    }
}