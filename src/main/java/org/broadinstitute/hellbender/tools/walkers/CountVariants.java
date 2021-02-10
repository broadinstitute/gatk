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
 * Count variant records in a VCF file, regardless of filter status. The tool gives the count at end of the standard out
 * (and optionally write it in a file).
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
    static final String USAGE_SUMMARY = "T Count variant records in a VCF file, regardless of filter status. The tool gives the count at end of the standard out\n" +
            " (and optionally write it in a file).\n" +
            "\n" +
            " <h3> Input </h3>\n" +
            " <ul>\n" +
            "     <li> A single VCF file. </li>\n" +
            " </ul>\n" +
            "\n" +
            " <h3> Usage example: </h3>\n" +
            " <pre>\n" +
            "  gatk CountVariants \\\n" +
            "      -V input_variants.vcf\n" +
            " </pre>";

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
