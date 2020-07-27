package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

import java.io.File;

@CommandLineProgramProperties(summary = "Compare Reference and Array VCFs for imputation pipeline. Adds variants either if not present" +
        "in the reference panel or if present in the reference panel and the reference and alternate alleles match OR if there is a simple" +
        "allele swap: panel's alternate allele matches array reference allele and panel's reference allele matches array variant allele. " +
        "Any multiallelic (in either reference panel or array) will be filtered out.",
        oneLineSummary = "Compare Reference and Array VCFs for imputation pipeline", programGroup = VariantManipulationProgramGroup.class)
public class MatchWithReferencePanel extends VariantWalker {

    @Argument(fullName="reference_panel",
            shortName = "panel", doc="Reference panel to match reference and alternate sites with",
            optional=false)
    private FeatureInput<VariantContext> referencePanel = null;


    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output vcf", optional=false)
    private File outputVcf = null;


    private VariantContextWriter vcfWriter;



    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(inputHeader);

        super.onTraversalStart();
    }


        @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {

        // Don't add array variant if it is multiallelic
        if (! vc.isBiallelic()) {
            logger.warn("More than 1 alternate allele in array variant for " + vc.getContig() + ":" + vc.getStart() +
                    ". The alternate alleles are " + vc.getAlternateAlleles().toString());
        } else {
            // add array variant it if isn't in the reference panel
            if (featureContext.getValues(referencePanel).size() == 0) {
                logger.warn("Reference panel does not contain variant at " + vc.getContig() + ":" + vc.getStart());
                vcfWriter.add(vc);
            } else {
                // panel should only contain biallelic sites
                if (featureContext.getValues(referencePanel).get(0).getAlternateAlleles().size() > 1) {
                    logger.warn("More than 1 alternate allele in reference panel for " + vc.getContig() + ":" + vc.getStart() +
                            ". The alternate alleles are " + featureContext.getValues(referencePanel).get(0).getAlternateAlleles().toString());
                } else {

                    Allele panelRef = featureContext.getValues(referencePanel).get(0).getReference();
                    Allele panelAlt = featureContext.getValues(referencePanel).get(0).getAlternateAllele(0);
                    Allele variantRef = vc.getReference();
                    Allele variantAlt = vc.getAlternateAllele(0);

                    // add array variant if it matches reference panel or is a simple swap between alternate and reference alleles
                    if ((panelRef == variantRef && panelAlt == variantAlt) || (panelRef.basesMatch(variantAlt) && panelAlt.basesMatch(variantRef))) {
                        vcfWriter.add(vc);
                    }
                }
            }
        }
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
