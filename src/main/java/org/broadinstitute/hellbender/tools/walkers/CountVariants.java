package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

/**
 *
 * Count variant records in a VCF file, regardless of filter status. The tool gives the count at end of the standard out.
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> A single VCF file. </li>
 * </ul>
 *
 * <h3> Usage example: </h3>
 * <pre>
 *  gatk CountVariants \
 *      -V input_variants.vcf
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = CountVariants.USAGE_SUMMARY,
        oneLineSummary = CountVariants.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)
public final class CountVariants extends VariantWalker{
    private long count = 0;

    static final String USAGE_ONE_LINE_SUMMARY = "Counts variant records in a VCF file, regardless of filter status.";
    static final String USAGE_SUMMARY = "This tool counts the variant records in a VCF file, regardless of filter status. " +
            "Because it counts the number of rows in the VCF, it does not necessarily reflect the number of variant " +
            "alleles. The count is returned at the end of the standard out.";

    @Override
    public void apply( VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
        count++;
    }

    @Override
    public Object onTraversalSuccess() {
        return count;
    }
}
