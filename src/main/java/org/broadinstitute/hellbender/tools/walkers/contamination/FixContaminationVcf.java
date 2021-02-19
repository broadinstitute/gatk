package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;


@CommandLineProgramProperties(
        summary = "Fix",
        oneLineSummary = "this now",
        programGroup = VariantFilteringProgramGroup.class
)
public class FixContaminationVcf extends VariantWalker { // or, the bed
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName =StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="", optional=false)
    private final String outputVcf = null;

    @Argument(fullName = "output2",
            doc="", optional=false)
    private final String badSitesVcf = null;

    private VariantContextWriter vcfWriter;
    private VariantContextWriter badSitesWriter;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart(){
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(getHeaderForVariants());

        badSitesWriter = createVCFWriter(new File(badSitesVcf));
        badSitesWriter.writeHeader(getHeaderForVariants());
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (! variant.getReference().basesMatch(referenceContext.getBases())){
            badSitesWriter.add(variant);
        } else {
            vcfWriter.add(variant);
        }

    }

    @Override
    public Object onTraversalSuccess(){
        vcfWriter.close();
        badSitesWriter.close();
        return "YEAAA BABY";
    }

}
