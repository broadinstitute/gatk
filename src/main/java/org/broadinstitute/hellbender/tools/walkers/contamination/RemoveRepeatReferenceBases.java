package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReferenceBases;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
@CommandLineProgramProperties(
        summary = "ah",
        oneLineSummary = "ah",
        programGroup = VariantFilteringProgramGroup.class
)
public class RemoveRepeatReferenceBases extends VariantWalker {
    private VariantContextWriter vcfWriter;
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName =StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file", optional=false)
    private final String outputVcf = null;

    @Override
    public void onTraversalStart(){
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(getHeaderForVariants());
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final String referenceBases = referenceContext.getKmerAround(variant.getStart(), 1);
        if (!(referenceBases.contains("AA") || referenceBases.contains("CC") ||
                referenceBases.contains("GG") || referenceBases.contains("TT"))){
            vcfWriter.add(variant);
        }
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}
