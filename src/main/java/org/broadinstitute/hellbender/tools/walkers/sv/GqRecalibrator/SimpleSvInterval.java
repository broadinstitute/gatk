package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

public class SimpleSvInterval implements Locatable {
    final SimpleInterval simpleInterval;
    final String svType;
    final int svLen;
    static final int BREAK_END_HALF_WIDTH = 121;
    static final Set<String> BREAK_END_SV_TYPES = new HashSet<>(Arrays.asList("BND", "CTX"));
    private static final Set<String> LARGE_CNV_TYPES = new HashSet<>(Arrays.asList("CNV", "DEL", "DUP"));
    static final int LARGE_CNV_MIN_SIZE = 5000;
    static final String SV_TYPE_KEY = "SVTYPE";
    static final String SV_LEN_KEY = "SVLEN";
    static final String CHR2_KEY = "CHR2";
    static final String CPX_INTERVALS_KEY = "CPX_INTERVALS";

    SimpleSvInterval(final SimpleInterval simpleInterval, final String svType, final int svLen) {
        this.simpleInterval = simpleInterval;
        this.svType = svType;
        this.svLen = svLen;
    }
    SimpleSvInterval(final String contig, final int start, final int end, final String svType, final int svLen) {
        this(new SimpleInterval(contig, start, end), svType, svLen);
    }

    static Stream<SimpleSvInterval> streamFrom(final VariantContext variantContext) {
        final String svType = variantContext.getAttributeAsString(SV_TYPE_KEY, null);
        if(svType == null) {
            throw new IllegalArgumentException(SV_TYPE_KEY + " must be defined");
        }
        if(BREAK_END_SV_TYPES.contains(svType)) {
            // Create two zero-length intervals
            final String chr2 = variantContext.getAttributeAsString(CHR2_KEY, variantContext.getContig());
            return Stream.of(
                new SimpleSvInterval(variantContext.getContig(), variantContext.getStart(), variantContext.getStart(),
                                     svType, 0),
                new SimpleSvInterval(chr2, variantContext.getEnd(), variantContext.getEnd(),
                                     svType, 0)
            );
        } else {
            Stream.Builder<SimpleSvInterval> simpleSvIntervalStreamBuilder = Stream.builder();

            final int svLen = variantContext.getAttributeAsInt(SV_LEN_KEY, Integer.MIN_VALUE);
            // Always add the primary interval
            simpleSvIntervalStreamBuilder.add(
                    new SimpleSvInterval(variantContext.getContig(), variantContext.getStart(), variantContext.getEnd(),
                                         svType, svLen)
            );
            // If there are CPX_INTERVALS, add them
            final List<String> cpxIntervals = variantContext.getAttributeAsStringList(CPX_INTERVALS_KEY, null);
            if(cpxIntervals != null) {
                cpxIntervals.forEach(cpxIntervalStr -> simpleSvIntervalStreamBuilder.add(getDistalTarget(cpxIntervalStr)));
            }
            return simpleSvIntervalStreamBuilder.build();
        }
    }

    private static SimpleSvInterval getDistalTarget(final String targetString) {
        final int underscoreIndex = targetString.indexOf("_");
        final SimpleInterval simpleInterval = new SimpleInterval(targetString.substring(underscoreIndex + 1));
        final String svType = targetString.substring(0, underscoreIndex);
        final int svLen = simpleInterval.getEnd() - simpleInterval.getStart();
        return new SimpleSvInterval(simpleInterval, svType, svLen);
    }

    static int getOrInferSvLen(final VariantContext variantContext) {
        final int svLen = variantContext.getAttributeAsInt(SV_LEN_KEY, Integer.MIN_VALUE);
        if(svLen == Integer.MIN_VALUE) {
            return streamFrom(variantContext).mapToInt(SimpleSvInterval::getOrInferSvLen).sum();
        } else {
            return svLen;
        }
    }

    @Override
    public String getContig() {
        return simpleInterval.getContig();
    }
    @Override
    public int getStart() {
        return simpleInterval.getStart();
    }
    @Override
    public int getEnd() {
        return simpleInterval.getEnd();
    }

    public int getOrInferSvLen() {
        if(svLen == Integer.MIN_VALUE) {
            if(BREAK_END_SV_TYPES.contains(svType)) {
                return 0;
            } else {
                return Integer.max(0, getEnd() - getStart());
            }
        } else {
            return svLen;
        }
    }

    Stream<SimpleInterval> streamGenomeTrackOverlapLocations() {
        if(LARGE_CNV_TYPES.contains(svType) && getOrInferSvLen() >= LARGE_CNV_MIN_SIZE) {
            // For large CNVs, check the primary interval for genome track overlaps
            return Stream.of(this.simpleInterval);
        } else {
            // For everything else, check the breakpoints
            return streamExpandedBreakEndIntervals();
        }
    }

    Stream<SimpleInterval> streamExpandedBreakEndIntervals() {
        if(getEnd() - getStart() <= 2 * BREAK_END_HALF_WIDTH) {
            return Stream.of(
                new SimpleInterval(getContig(), getStart() - BREAK_END_HALF_WIDTH, getEnd() + BREAK_END_HALF_WIDTH)
            );
        } else {
            return Stream.of(
                new SimpleInterval(getContig(), getStart() - BREAK_END_HALF_WIDTH, getStart() + BREAK_END_HALF_WIDTH),
                new SimpleInterval(getContig(), getEnd() - BREAK_END_HALF_WIDTH, getEnd() + BREAK_END_HALF_WIDTH)
            );
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SimpleSvInterval that = (SimpleSvInterval) o;
        return this.simpleInterval.equals(that.simpleInterval) &&
               this.svType.equals(that.svType) &&
               this.svLen == that.svLen;
    }

    @Override
    public int hashCode() {
        return 31 * (this.simpleInterval.hashCode() + this.svType.hashCode() + this.svLen);
    }
}
