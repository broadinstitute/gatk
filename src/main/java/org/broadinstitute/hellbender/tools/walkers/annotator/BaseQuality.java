package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.List;
import java.util.OptionalDouble;
import java.util.OptionalInt;

/**
 * Median base quality of bases supporting each allele.
 *
 * Created by David Benjamin on 3/20/17.
 */
public class BaseQuality extends PerAlleleAnnotation {
    public static final String KEY = "MBQ";

    @Override
    protected int aggregate(final List<Integer> values) {
        return values.isEmpty() ? 0 : GATKProtectedMathUtils.median(Ints.toArray(values));
    }

    @Override
    protected String getVcfKey() { return KEY; }

    @Override
    protected String getDescription() { return "median base quality"; }

    @Override
    protected boolean includeRefAllele() { return true; }

    @Override
    protected OptionalInt getValueForRead(final GATKRead read, final VariantContext vc) {
        Utils.nonNull(read);

        final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), vc.getStart(), ReadUtils.ClippingTail.RIGHT_TAIL, true);
        return offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED || AlignmentUtils.isInsideDeletion(read.getCigar(), offset) ?
                OptionalInt.empty() : OptionalInt.of(read.getBaseQuality(offset));
    }
}
