package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;
import java.util.Set;

/**
 * A {@link VariantContextWriter} decorator which filters out variants that don't match a given set of intervals.
 */
public class IntervalFilteringVcfWriter implements VariantContextWriter {

    /**
     * Comparison modes which allow matching intervals in different ways.
     */
    public enum Mode {

        /**
         * Matches if the query starts within any of the given intervals.
         */
        STARTS_IN{
            @Override
            boolean test(OverlapDetector<? extends Locatable> detector, final VariantContext query) {
                final SimpleInterval startPosition = new SimpleInterval(query.getContig(), query.getStart(), query.getStart());
                return detector.overlapsAny(startPosition);
            }
        },

        /**
         * Matches if the query ends within any of the given intervals
         */
        ENDS_IN{
            @Override
            boolean test(final OverlapDetector<? extends Locatable> detector, final VariantContext query) {
                final SimpleInterval endPosition = new SimpleInterval(query.getContig(), query.getEnd(), query.getEnd());
                return detector.overlapsAny(endPosition);
            }
        },

        /**
         * Matches if any part of the query overlaps any one of the given intervals
         */
        OVERLAPS{
            @Override
            boolean test(final OverlapDetector<? extends Locatable> detector, final VariantContext query) {
                return detector.overlapsAny(query);
            }
        },

        /**
         * Matches if the entirety of the query is contained within one of the intervals
         */
        CONTAINED {
            @Override
            boolean test(final OverlapDetector<? extends Locatable> detector, final VariantContext query) {
                final Set<? extends Locatable> overlaps = detector.getOverlaps(query);
                for( final Locatable loc : overlaps){
                    if(loc.contains(query)){
                        return true;
                    }
                }
                return false;
            }
        },

        /**
         * Always matches, may be used to not perform any filtering, alternatively a
         */
        ANYWHERE {
            @Override
            boolean test(final OverlapDetector<? extends Locatable> detector, final VariantContext query) {
                return true;
            }
        };

        /**
         * @param detector The OverlapDetector to compare against
         * @param query The variant being tested
         * @return true iff the variant matches the given intervals
         */
        abstract boolean test(OverlapDetector<? extends Locatable> detector, VariantContext query);
    }

    private final VariantContextWriter writer;
    private final OverlapDetector<SimpleInterval> detector;
    private final Mode mode;

    /**
     * @param writer the writer to wrap
     * @param intervals the intervals to compare against, note that these are not merged so if they should be merged than the input list should be preprocessed
     * @param mode the matching mode to use
     */
    public IntervalFilteringVcfWriter(final VariantContextWriter writer, List<SimpleInterval> intervals, Mode mode) {
        Utils.nonNull(writer);
        Utils.nonEmpty(intervals);
        Utils.nonNull(mode);

        this.writer = writer;
        this.detector = OverlapDetector.create(intervals);
        this.mode = mode;
    }

    @Override
    public void writeHeader(final VCFHeader header) {
        writer.writeHeader(header);
    }

    @Override
    public void setHeader(final VCFHeader header) {
        writer.setHeader(header);
    }

    @Override
    public void close() {
        writer.close();
    }

    @Override
    public boolean checkError() {
        return writer.checkError();
    }

    /**
     * Add the given variant to the writer and output it if it matches.
     * @param vc the variant to potentially write
     */
    @Override
    public void add(final VariantContext vc) {
        if(mode.test(detector, vc)) {
            writer.add(vc);
        }
    }

}
