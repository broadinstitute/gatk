package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.TwoPassVariantWalker;
import org.ojalgo.array.blas.COPY;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.QUAL_BY_DEPTH_KEY;

@CommandLineProgramProperties(
        summary = "Example variant walker that makes two passes through a vcf",
        oneLineSummary = "Example two variant walker",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public class ExampleTwoPassVariantWalker extends TwoPassVariantWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output vcf", optional = true)
    private File outputVcf = null;

    final List<Double> qualByDepths = new ArrayList<>();

    private VariantContextWriter vcfWriter;

    public static String QD_P_VALUE_KEY_NAME = "P_VALUE";

    public static String COPY_OF_QD_KEY_NAME = "QD_COPY";

    private double averageQualByDepth;

    private double varianceOfQDs;

    int counter = 0;

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        inputHeader.addMetaDataLine(new VCFInfoHeaderLine(QD_P_VALUE_KEY_NAME, 1, VCFHeaderLineType.Float, "p value of the quality by depth"));
        inputHeader.addMetaDataLine(new VCFInfoHeaderLine(COPY_OF_QD_KEY_NAME, 1, VCFHeaderLineType.Float, "copy of the QD INFO field"));
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(inputHeader);
    }

    @Override
    protected void firstPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        qualByDepths.add(variant.getAttributeAsDouble(QUAL_BY_DEPTH_KEY, 0.0));
    }

    @Override
    protected void secondPassApply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        VariantContextBuilder vcb = new VariantContextBuilder(variant);
        final double qualByDepth = qualByDepths.get(counter++);
        final double pValue = 0.05;
        vcb.attribute(QD_P_VALUE_KEY_NAME, pValue);
        vcb.attribute(COPY_OF_QD_KEY_NAME, qualByDepth);

        vcfWriter.add(vcb.make());
    }

    @Override
    protected void afterFirstPass() {
        final int n = qualByDepths.size();
        if (n == 1){
            return;
        }
        averageQualByDepth = qualByDepths.stream().mapToDouble(x -> x).average().orElse(0.0);
        varianceOfQDs = (1/(n-1))* qualByDepths.stream().mapToDouble(x -> Math.pow((x - averageQualByDepth), 2.0)).sum();
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
