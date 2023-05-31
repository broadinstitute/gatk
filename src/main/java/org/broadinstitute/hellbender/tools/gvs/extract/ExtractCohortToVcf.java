package org.broadinstitute.hellbender.tools.gvs.extract;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.utils.io.IOUtils;


@SuppressWarnings("unused")
@CommandLineProgramProperties(
        summary = "(\"ExtractCohortToVcf\") - Filter and extract variants from BigQuery to a VCF output.",
        oneLineSummary = "Tool to extract variants from BigQuery to a VCF output for a subset of samples.",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ExtractCohortToVcf extends ExtractCohort {
    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which annotated variants should be written."
    )
    protected String outputVcfPathString = null;

    protected VariantContextWriter vcfWriter = null;

    @Override
    protected void onStartup() {
        super.onStartup();

        vcfWriter = createVCFWriter(IOUtils.getPath(outputVcfPathString));
        vcfWriter.writeHeader(header);
    }

    @Override
    protected void apply(VariantContext variantContext) {
        if (variantContext != null) {
            // Add the variant contexts that aren't filtered or add everything if we aren't excluding anything
            if (variantContext.isNotFiltered() || !excludeFilteredSites) {
                vcfWriter.add(variantContext);
            }
            progressMeter.update(variantContext);
        }
    }

    @Override
    protected void onShutdown() {
        super.onShutdown();

        // Close up our writer if we have to:
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
