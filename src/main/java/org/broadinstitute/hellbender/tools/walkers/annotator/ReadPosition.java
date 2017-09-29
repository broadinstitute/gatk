package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.List;
import java.util.OptionalInt;

/**
 * Median distance of variant starts from ends of reads supporting each allele.
 *
 * Created by David Benjamin on 3/20/17.
 */
public class ReadPosition extends PerAlleleAnnotation {
    public static final String KEY = "MPOS";

    @Override
    protected int aggregate(final List<Integer> values) {
        return values.isEmpty() ? 0 : GATKProtectedMathUtils.median(Ints.toArray(values));
    }

    @Override
    protected String getVcfKey() { return KEY; }

    @Override
    protected String getDescription() { return "median distance from end of read"; }

    @Override
    protected OptionalInt getValueForRead(final GATKRead read, final VariantContext vc) {
        Utils.nonNull(read);
        final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), vc.getStart(), ReadUtils.ClippingTail.RIGHT_TAIL, true);
        if ( offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED || AlignmentUtils.isInsideDeletion(read.getCigar(), offset)) {
            return OptionalInt.empty();
        }

        int readPosition = AlignmentUtils.calcAlignmentByteArrayOffset(read.getCigar(), offset, false, 0, 0);
        final int numAlignedBases = AlignmentUtils.getNumAlignedBasesCountingSoftClips( read );
        final int distanceFromEnd = Math.min(readPosition, numAlignedBases - readPosition - 1);
        return OptionalInt.of(distanceFromEnd);
    }
}