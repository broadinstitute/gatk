package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

/**
 *
 * Count variants in a VCF file.
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> A single VCF file. </li>
 * </ul>
 *
 * <h3> Example: </h3>
 * <pre>
 *  gatk CountVariants \
 *      -V input_variants.vcf
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = CountVariants.USAGE_SUMMARY,
        oneLineSummary = CountVariants.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantProgramGroup.class
)
public final class CountVariants extends VariantWalker{
    private long count = 0;

    static final String USAGE_ONE_LINE_SUMMARY = "Count variants in a VCF file";
    static final String USAGE_SUMMARY = "Walks over the input data set, calculating the number of variants seen.";

    @Override
    public void apply( VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
        count++;
    }

    @Override
    public Object onTraversalSuccess() {
        return count;
    }
}
