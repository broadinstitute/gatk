package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "Functional Annotation Filtration",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class FilterFuncotations extends VariantWalker {
    private static Logger logger = LogManager.getLogger(FilterFuncotations.class);
    //==================================================================================================================
    // Arguments:

    //-----------------------------------------------------
    // Required args:

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which annotated variants should be written.")
    protected File outputFile;


    private VariantContextWriter outputVcfWriter;
    private String[] funcotationKeys;
    private List<FuncotationFilter> funcotationFilters = new ArrayList<>();

    @Override
    public void onTraversalStart() {
        registerFilters();
        VCFHeader vcfHeader = getHeaderForVariants();
        final VCFInfoHeaderLine funcotationHeaderLine = vcfHeader.getInfoHeaderLine(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME);
        if (funcotationHeaderLine != null) {
            funcotationKeys = FuncotatorUtils.extractFuncotatorKeysFromHeaderDescription(funcotationHeaderLine.getDescription());
            outputVcfWriter = createVCFWriter(outputFile);
            funcotationFilters.forEach(filter -> vcfHeader.addMetaDataLine(new VCFFilterHeaderLine(filter.getFilterName())));
            outputVcfWriter.writeHeader(vcfHeader);
        } else {
            logger.error("Input VCF does not have Funcotator annotations.");
        }
    }

    private void registerFilters() {
        funcotationFilters.add(new ClinVarFilter());
        funcotationFilters.add(new LofFilter());
        funcotationFilters.add(new LmmFilter());
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        outputVcfWriter.add(applyFilters(variant));
    }

    private VariantContext applyFilters(VariantContext variant) {
        Map<Allele, FuncotationMap> funcs = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                funcotationKeys, variant, "test", "FilterFuncs"
        );
        VariantContextBuilder variantContextBuilder = new VariantContextBuilder(variant);
        funcs.values().forEach(funcotationMap ->
                funcotationMap.getTranscriptList().forEach(transcriptId -> {
                    funcotationFilters.forEach(filter -> {
                        if (filter.checkFilter(funcotationMap.get(transcriptId))) {
                            variantContextBuilder.filter(filter.getFilterName());
                        }
                    });
                }));
        return variantContextBuilder.make();
    }


    @Override
    public void closeTool() {
        if (outputVcfWriter != null) {
            outputVcfWriter.close();
        }
    }
}

class ClinVarFilter extends FuncotationFilter {

    ClinVarFilter() {
        super("CLINVAR");
    }

    @Override
    public Boolean checkFilter(List<Funcotation> funcotations) {
        return false;
    }
}

class LofFilter extends FuncotationFilter {

    LofFilter() {
        super("LOF");
    }

    @Override
    public Boolean checkFilter(List<Funcotation> funcotations) {
        return false;
    }
}

class LmmFilter extends FuncotationFilter {
    LmmFilter() {
        super("LMM");
    }


    @Override
    public Boolean checkFilter(List<Funcotation> funcotations) {
        return false;
    }
}

abstract class FuncotationFilter {

    private final String filterName;

    FuncotationFilter(String filterName) {
        this.filterName = filterName;
    }

    public abstract Boolean checkFilter(List<Funcotation> funcotations);

    public String getFilterName() {
        return filterName;
    }
}