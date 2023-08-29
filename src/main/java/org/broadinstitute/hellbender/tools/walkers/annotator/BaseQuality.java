package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;
import java.util.OptionalDouble;
import java.util.OptionalInt;

/**
 * Median base quality of bases supporting each allele.
 *
 * <p>The output is an array containing, for each allele, the median base quality at the variant (for SNVs) and one base before the variant (for indels) over all reads that best match that allele.</p>
 *
 * <p>For example, for variant context with ref allele A and alt allele C we use base qualities at the A and C.  For variant context with ref allele AG and alt allele A (deletion),
 * we use base qualities at the A.  For variant context with ref allele A and alt allele AG (insertion) we use base qualities at the A.</p>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Median base quality of bases supporting each allele (MBQ)")
public class BaseQuality extends PerAlleleAnnotation implements StandardMutectAnnotation {

    @Override
    protected int aggregate(final List<Integer> values) {
        return values.isEmpty() ? 0 : MathUtils.median(Ints.toArray(values));
    }

    @Override
    protected String getVcfKey() { return GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY; }

    @Override
    protected boolean includeRefAllele() { return true; }

    @Override
    protected OptionalInt getValueForRead(final GATKRead read, final VariantContext vc) {
        return getBaseQuality(read, vc);
    }

    public static OptionalInt getBaseQuality(final GATKRead read, final VariantContext vc) {
        if (vc.getStart() < read.getStart() || read.getEnd() < vc.getStart()) {
            return OptionalInt.empty();
        }
        final OptionalDouble result = BaseQualityRankSumTest.getReadBaseQuality(read, vc);
        return result.isPresent() ? OptionalInt.of((int) FastMath.round(result.getAsDouble())) : OptionalInt.empty();
    }
}
