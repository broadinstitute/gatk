package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

public final class SVStratification {

    final GATKSVVCFConstants.StructuralVariantAnnotationType svType;
    final int minSize;  // inclusive
    final int maxSize;  // exclusive
    final String contextName;
    final OverlapDetector<Locatable> contextIntervals;
    final String name;

    public SVStratification(final GATKSVVCFConstants.StructuralVariantAnnotationType svType,
                            final int minSize, final int maxSize, final String contextName,
                            final List<Locatable> contextIntervals) {
        Utils.nonNull(contextName);
        Utils.nonNull(contextIntervals);
        if (maxSize >= 0 && minSize >= 0 && maxSize < minSize) {
            throw new IllegalArgumentException("Max size cannot be less than min size");
        }
        this.svType = svType;
        // Map min from any negative number to -1
        final String minSizeString;
        if (minSize < 0) {
            this.minSize = -1;
            minSizeString = "";
        } else {
            this.minSize = minSize;
            minSizeString = "_" + minSize;
        }
        // Map max from any negative number to infinity
        final String maxSizeString;
        if (maxSize < 0) {
            this.maxSize = Integer.MAX_VALUE;
            maxSizeString = "";
        } else {
            this.maxSize = maxSize;
            maxSizeString = "_" + maxSize;
        }
        this.contextName = contextName;
        this.contextIntervals = OverlapDetector.create(contextIntervals);
        this.name = svType.name() + minSizeString + maxSizeString + "_" + contextName;
    }

    public boolean matches(final VariantContext variant, final double overlapFraction,
                           final int numBreakpointOverlaps, final SAMSequenceDictionary dictionary) {
        final SVCallRecord record = SVCallRecordUtils.create(variant, dictionary);
        return matchesType(record) && matchesSize(record) && matchesContext(record, overlapFraction, numBreakpointOverlaps);
    }

    protected boolean matchesType(final SVCallRecord record) {
        return record.getType() == svType;
    }

    protected boolean matchesSize(final SVCallRecord record) {
        final Integer length = record.getLength();
        return length == null || (length >= minSize && length < maxSize);
    }

    public boolean matchesContext(final SVCallRecord record,
                                  final double overlapFraction,
                                  final int numBreakpointOverlaps) {
        Utils.nonNull(record);
        Utils.validate(overlapFraction >= 0 && overlapFraction <= 1,
                "Invalid overlap fraction threshold " + overlapFraction);
        Utils.validate(numBreakpointOverlaps >= 0 && numBreakpointOverlaps <= 2,
                "Invalid breakpoint overlaps threshold " + numBreakpointOverlaps);
        Utils.validate(!(overlapFraction == 0 && numBreakpointOverlaps == 0),
                "Overlap fraction and overlapping breakpoints thresholds cannot both be 0");
        if (record.isIntrachromosomal()) {
            return matchesContextIntrachromosomal(record, overlapFraction, numBreakpointOverlaps);
        } else {
            return matchesContextBreakpointOverlap(record, numBreakpointOverlaps);
        }
    }

    protected boolean matchesContextIntrachromosomal(final SVCallRecord record,
                                                     final double overlapFraction,
                                                     final int numBreakpointOverlaps) {
        return matchesContextOverlapFraction(record, overlapFraction) && matchesContextBreakpointOverlap(record, numBreakpointOverlaps);
    }

    protected boolean matchesContextOverlapFraction(final SVCallRecord record, final double overlapFraction) {
        // TODO handle insertions / cpx
        if (overlapFraction > 0) {
            final SimpleInterval interval = new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionB());
            long overlapLength = 0;
            for (final Locatable overlap : contextIntervals.getOverlaps(interval)) {
                overlapLength += interval.intersect(overlap).size();
            }
            return overlapLength / (double) interval.getLengthOnReference() >= overlapFraction;
        } else {
            return true;
        }
    }

    protected boolean matchesContextBreakpointOverlap(final SVCallRecord record, final int numBreakpointOverlaps) {
        if (numBreakpointOverlaps > 0) {
            final SimpleInterval intervalA = new SimpleInterval(record.getContigA(), record.getPositionA(), record.getPositionA());
            final SimpleInterval intervalB = new SimpleInterval(record.getContigB(), record.getPositionB(), record.getPositionB());
            return countAnyContextOverlap(intervalA) + countAnyContextOverlap(intervalB) >= numBreakpointOverlaps;
        } else {
            return true;
        }
    }

    protected int countAnyContextOverlap(final SimpleInterval interval) {
        if (contextIntervals.overlapsAny(interval)) {
            return 1;
        } else {
            return 0;
        }
    }

    public boolean isMutuallyExclusive(final SVStratification other) {
        Utils.nonNull(other);
        if (svType != other.getSvType()) {
            return true;
        } else if (!contextName.equals(other.getContextName())) {
            return true;
        } else if (minSize >= other.getMaxSize() || other.getMinSize() >= maxSize) {
            return true;
        } else {
            return false;
        }
    }

    public GATKSVVCFConstants.StructuralVariantAnnotationType getSvType() {
        return svType;
    }

    public Integer getMinSize() {
        return minSize;
    }

    public Integer getMaxSize() {
        return maxSize;
    }

    public String getContextName() {
        return contextName;
    }

    public String getName() {
        return name;
    }

    @Override
    public String toString() {
        return "SVStratification{" +
                "svType=" + svType +
                ", minSize=" + minSize +
                ", maxSize=" + maxSize +
                ", contextName='" + contextName + '\'' +
                ", name='" + name + '\'' +
                '}';
    }
}
