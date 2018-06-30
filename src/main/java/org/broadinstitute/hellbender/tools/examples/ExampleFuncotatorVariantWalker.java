package org.broadinstitute.hellbender.tools.examples;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;

import java.util.Map;

/**
 * TODO: Docs
 * TODO: Integration tests
 */
@CommandLineProgramProperties(
        summary = "Example variant walker that reads a funcotated VCF and prints out a summary.",
        oneLineSummary = "Example variant walker for a funcotated VCF",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public class ExampleFuncotatorVariantWalker extends VariantWalker {
    private String[] funcotationFieldNames;

    @Override
    public void onTraversalStart() {
        funcotationFieldNames = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(getHeaderForVariants()
            .getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME).getDescription());
        logger.info(StringUtils.join(funcotationFieldNames, ","));
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final Map<Allele, FuncotationMap> alleleMap =  FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(funcotationFieldNames, variant, "Gencode_27_annotationTranscript", "FAKE");
    }
}
