package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalTextOutputArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

/**
 *
 * Count variant records in a VCF file, regardless of filter status. The tool prints the count to standard output
 * (and can optionally write it to a file).
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
public final class CountVariants extends VariantWalker {
    private long count = 0;

    static final String USAGE_ONE_LINE_SUMMARY = "Counts variant records in a VCF file, regardless of filter status.";
    static final String USAGE_SUMMARY = "This tool counts the variant records in a VCF file, regardless of filter status. " +
            "Because it counts the number of rows in the VCF, it does not necessarily reflect the number of variant " +
            "alleles. The count is printed to standard output (and may optionally be written to a file)";

    @ArgumentCollection
    final public OptionalTextOutputArgumentCollection out = new OptionalTextOutputArgumentCollection();

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        count++;
    }

    @Override
    public Object onTraversalSuccess() {
        out.print(count);
        return count;
    }
}
