package org.broadinstitute.hellbender.tools.exome;


import htsjdk.tribble.bed.BEDFeature;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Utilities to pad targets, so tadded targets are not allowed to overlap.  Basically, if padded targets would overlap,
 * the padding goes to the midpoint between two targets.
 */
public final class TargetPadder {

    private TargetPadder() {}

    /**
     * Create a TargetCollection of padded Targets from BEDFeatures.
     *
     * @param inputTargets BEDFeatures for padding.  Cannot be {@code null}.
     * @param paddingInBases Number of bases to pad.  Must be greater than or equal to zero.
     * @return never {@code null}
     */
    public static TargetCollection<Target> padTargetsFromBEDFeatures(final TargetCollection<? extends BEDFeature> inputTargets, final int paddingInBases) {
        ParamUtils.isPositiveOrZero(paddingInBases, "Padding must be >= 0");
        Utils.nonNull(inputTargets, "Input targets cannot be null.");
        final List<Target> targetList = inputTargets.targets().stream().map(t -> createTargetFromBEDFeature(t)).collect(Collectors.toList());
        return padTargets(new HashedListTargetCollection<>(targetList), paddingInBases);
    }
    /**
     * Create a TargetCollection of padded Targets from Targets.
     *
     * Note: Do not use this method on subclasses such as TargetCoverage.
     *
     * @param inputTargets BEDFeatures for padding.  Cannot be {@code null}.
     * @param paddingInBases Number of bases to pad.  Must be greater than or equal to zero.
     * @return never {@code null}
     */
    public static TargetCollection<Target> padTargets(final TargetCollection<Target> inputTargets, final int paddingInBases) {
        ParamUtils.isPositiveOrZero(paddingInBases, "Padding must be >= 0");
        Utils.nonNull(inputTargets, "Input targets cannot be null.");
        final List<Target> paddedTargets = inputTargets.targets().stream().map(t -> createPaddedTarget(t, paddingInBases)).collect(Collectors.toList());

        // Alter the padded Targets to eliminate overlaps
        for (int i = 0; i < paddedTargets.size() - 1; i++){
            final Target thisTarget = paddedTargets.get(i);
            final Target nextTarget = paddedTargets.get(i + 1);
            if (thisTarget.getInterval().overlaps(nextTarget.getInterval())) {

                final int originalThisEnd = inputTargets.targets().get(i).getEnd();
                final int originalNextStart = inputTargets.targets().get(i + 1).getStart();

                final int newThisEnd = (originalThisEnd + originalNextStart)/2;
                final int newNextStart = newThisEnd + 1;

                paddedTargets.set(i, new Target(thisTarget.getName(), new SimpleInterval(thisTarget.getContig(), thisTarget.getInterval().getStart(), newThisEnd)));
                paddedTargets.set(i+1, new Target(nextTarget.getName(), new SimpleInterval(nextTarget.getContig(), newNextStart, nextTarget.getEnd())));
            }

        }

        return new HashedListTargetCollection<>(paddedTargets);
    }

    protected static Target createPaddedTarget(final Target t, final int paddingInBases){
        return new Target(t.getName(), new SimpleInterval(t.getContig(), Math.max(1, t.getStart() - paddingInBases), t.getEnd() + paddingInBases));
    }

    protected static Target createTargetFromBEDFeature(final BEDFeature bedFeature) {
        return new Target(bedFeature.getName(), new SimpleInterval(bedFeature.getContig(), bedFeature.getStart(), bedFeature.getEnd()));
    }
}
