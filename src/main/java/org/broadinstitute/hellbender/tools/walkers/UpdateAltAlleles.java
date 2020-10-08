package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.Level;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

@CommandLineProgramProperties(
        summary = "Tool that enables swapping alternate alleles that are coded as '.' with alternate alleles from another VCF",
        oneLineSummary = "Adds alternate alleles from another file",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public class UpdateAltAlleles extends VariantWalker {
  @Argument(fullName="referenceVCF", shortName="RP", doc="File from which to get alternate alleles (Reference Panel)")
    private FeatureInput<VariantContext> referenceVariants;


    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output vcf", optional = true)
    private File outputVcf = null;
    private VariantContextWriter vcfWriter;


    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        vcfWriter = createVCFWriter(outputVcf);
        vcfWriter.writeHeader(inputHeader);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        VariantContextBuilder vcb = new VariantContextBuilder(variant);

        if (! variant.isVariant()) {
            if (! featureContext.getValues(referenceVariants).isEmpty()) {
                for (VariantContext panel_variant:featureContext.getValues(referenceVariants)) {
                    if (panel_variant.getReference() == variant.getReference() && panel_variant.isSNP()) {
                        vcb.alleles(panel_variant.getAlleles());
                        vcfWriter.add(vcb.make());
                        break;
                    }
                }
            }
        } else {
            vcfWriter.add(vcb.make());
        }
    }


    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}

