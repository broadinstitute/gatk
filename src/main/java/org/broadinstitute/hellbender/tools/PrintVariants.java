package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

import java.util.Collections;

/**
 *  Example class to illustrate how to write a VariantWalker.
 */
@CommandLineProgramProperties(
        usage = "Prints variants from the input to the output.",
        usageShort = "Print variants",
        programGroup = VariantProgramGroup.class
)
public final class PrintVariants extends VariantWalker {

    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String out = null;

    private VariantContextWriter writer;

    @Override
    public void onTraversalStart() {
        writer = new VariantContextWriterBuilder().setOutputFile(out).setOutputFileType(VariantContextWriterBuilder.OutputType.VCF).unsetOption(Options.INDEX_ON_THE_FLY).build();
        writer.writeHeader(getHeaderForVariants());
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        writer.add(variant);
    }

    @Override
    public Object onTraversalDone() {
        CloserUtil.close(writer);
        return null;
    }
}
