package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.apache.commons.math3.stat.descriptive.rank.Median;


import java.util.List;
import java.util.OptionalDouble;
import java.util.OptionalInt;

/**
 * Median number of clipped (hard and high-quality soft) bases in reads supporting each allele.
 *
 * Created by David Benjamin on 3/20/17.
 */
public class ClippedBases extends PerAlleleAnnotation {
    private static byte BASE_QUAL_THRESHOLD = 25;

    public static final String KEY = "MCL";

    @Override
    protected int aggregate(final List<Integer> values) {
        return values.isEmpty() ? 0 : GATKProtectedMathUtils.median(Ints.toArray(values));
    }

    @Override
    protected String getVcfKey() { return KEY; }

    @Override
    protected String getDescription() { return "number of soft- and hard- clipped bases"; }

    @Override
    protected OptionalInt getValueForRead(final GATKRead read, final VariantContext vc) {
        Utils.nonNull(read);
        return OptionalInt.of(AlignmentUtils.getNumHardClippedBases(read) + AlignmentUtils.calcNumHighQualitySoftClips(read, BASE_QUAL_THRESHOLD));
    }
}