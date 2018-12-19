package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


/**
 * This walker makes two traversals through variants in a vcf. During the first pass, it collects the qual-by-depth (QD) annotation
 * of each variant and stores it in a list. After the first traversal, it computes the vcf-wide average and variance of QD's.
 * The walker will then traverse the variants the second time and annotate each variant with the distance of its QD from the mean
 * in units of standard deviations.
 *
 */
@CommandLineProgramProperties(
        summary = "Example variant walker that makes two passes through a vcf",
        oneLineSummary = "Example two-pass variant walker",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public class ExampleMultiplePassVariantWalker extends TwoPassVariantWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output vcf", optional = true)
    private File outputVcf = null;

    private final List<Double> qualByDepths = new ArrayList<>();

    private VariantContextWriter vcfWriter;

    public static final String QD_DISTANCE_FROM_MEAN = "QD_DIST";

    public static final String COPY_OF_QD_KEY_NAME = "QD_COPY";

    private double averageQualByDepth;

    private double sampleVarianceOfQDs;

    private int counter = 0;

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        inputHeader.addMetaDataLine(new VCFInfoHeaderLine(QD_DISTANCE_FROM_MEAN, 1, VCFHeaderLineType.Float, "distance from the average QD value in the units of standard deviations"));
        inputHeader.addMetaDataLine(new VCFInfoHeaderLine(COPY_OF_QD_KEY_NAME, 1, VCFHeaderLineType.Float, "copy of the QD INFO field"));
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(inputHeader);
    }

    @Override
    protected void firstPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        qualByDepths.add(variant.getAttributeAsDouble(GATKVCFConstants.QUAL_BY_DEPTH_KEY, 0.0));
    }

    @Override
    protected void afterFirstPass() {
        final int n = qualByDepths.size();
        if (n == 1){
            return;
        }
        averageQualByDepth = qualByDepths.stream().mapToDouble(x -> x).average().orElse(0.0);
        sampleVarianceOfQDs = (1.0/(n - 1.0))* qualByDepths.stream().mapToDouble(x -> Math.pow((x - averageQualByDepth), 2.0)).sum();
    }

    @Override
    protected void secondPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        final double qualByDepth = qualByDepths.get(counter++);
        final double distanceFromMean = Math.abs((qualByDepth - averageQualByDepth)/ Math.sqrt(sampleVarianceOfQDs));
        vcb.attribute(QD_DISTANCE_FROM_MEAN, distanceFromMean);
        vcb.attribute(COPY_OF_QD_KEY_NAME, qualByDepth);

        vcfWriter.add(vcb.make());
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
