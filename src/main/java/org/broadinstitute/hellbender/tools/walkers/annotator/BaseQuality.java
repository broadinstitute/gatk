package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.List;
import java.util.OptionalInt;

/**
 * Median base quality of bases supporting each alt allele.
 *
 * <p>The output is an array containing, for each alt allele, the median base quality at the variant (for SNVs) and one base before the variant (for indels) over all reads that best match that allele.</p>
 *
 * <p>For example, for variant context with ref allele A and alt allele C we use base qualities at the C.  For variant context with ref allele AG and alt allele A (deletion),
 * we use base qualities at the A.  For variant context with ref allele A and alt allele AG (insertion) we use base qualities at the A.</p>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Median base quality of bases supporting each allele (MBQ)")
public class BaseQuality extends PerAlleleAnnotation implements StandardMutectAnnotation {
    public static final String KEY = "MBQ";

    @Override
    protected int aggregate(final List<Integer> values) {
        return values.isEmpty() ? 0 : MathUtils.median(Ints.toArray(values));
    }

    @Override
    protected String getVcfKey() { return KEY; }

    @Override
    protected String getDescription() { return "median base quality"; }

    @Override
    protected OptionalInt getValueForRead(final GATKRead read, final VariantContext vc) {
        Utils.nonNull(read);

        final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), vc.getStart(), ReadUtils.ClippingTail.RIGHT_TAIL, true);
        return offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED || AlignmentUtils.isInsideDeletion(read.getCigar(), offset) ?
                OptionalInt.empty() : OptionalInt.of(read.getBaseQuality(offset));
    }
}
