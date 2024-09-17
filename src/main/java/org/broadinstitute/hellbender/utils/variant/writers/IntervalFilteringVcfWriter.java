package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineParser;
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
    public enum Mode implements CommandLineParser.ClpEnum {

        /**
         * Matches if the query starts within any of the given intervals.
         */
        STARTS_IN("starts within any of the given intervals"){

            @Override
            boolean test(final OverlapDetector<? extends Locatable> detector, final VariantContext query) {
                final SimpleInterval startPosition = new SimpleInterval(query.getContig(), query.getStart(), query.getStart());
                return detector.overlapsAny(startPosition);
            }
            @Override
            String getName() {return "STARTS_IN";}
        },

        /**
         * Matches if the query ends within any of the given intervals
         */
        ENDS_IN("ends within any of the given intervals"){
            @Override
            boolean test(final OverlapDetector<? extends Locatable> detector, final VariantContext query) {
                final SimpleInterval endPosition = new SimpleInterval(query.getContig(), query.getEnd(), query.getEnd());
                return detector.overlapsAny(endPosition);
            }
            @Override
            String getName() {return "ENDS_IN";}
        },

        /**
         * Matches if any part of the query overlaps any one of the given intervals
         */
        OVERLAPS("overlaps any of the given intervals"){
            @Override
            boolean test(final OverlapDetector<? extends Locatable> detector, final VariantContext query) {
                return detector.overlapsAny(query);
            }
            @Override
            String getName() {return "OVERLAPS";}
        },

        /**
         * Matches if the entirety of the query is contained within one of the intervals.  Note that adjacent intervals
         * may be merged into a single interval depending specified "--interval-merging-rule".
         */
        CONTAINED("contained completely within a contiguous block of intervals without overlap") {
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
            @Override
            String getName() {return "CONTAINED";}
        },

        /**
         * Always matches, may be used to not perform any filtering
         */
        ANYWHERE("no filtering") {
            @Override
            boolean test(final OverlapDetector<? extends Locatable> detector, final VariantContext query) {
                return true;
            }
            @Override
            String getName() {return "ANYWHERE";}
        };

        private final String doc;

        /**
         * @param detector The OverlapDetector to compare against
         * @param query The variant being tested
         * @return true iff the variant matches the given intervals
         */
        abstract boolean test(final OverlapDetector<? extends Locatable> detector, final VariantContext query);
        abstract String getName();

        private Mode(String doc){
            this.doc = doc;

        }

        @Override
        public String getHelpDoc() {
            return doc;
        }
    }

    private final VariantContextWriter underlyingWriter;
    private final OverlapDetector<SimpleInterval> detector;
    private final Mode mode;
    private static int filteredCount = 0;
    protected final Logger logger = LogManager.getLogger(this.getClass());

    /**
     * @param writer the writer to wrap
     * @param intervals the intervals to compare against, note that these are not merged so if they should be merged than the input list should be preprocessed
     * @param mode the matching mode to use
     */
    public IntervalFilteringVcfWriter(final VariantContextWriter writer, final List<SimpleInterval> intervals, final Mode mode) {
        Utils.nonNull(writer);
        Utils.nonEmpty(intervals);
        Utils.nonNull(mode);

        this.underlyingWriter = writer;
        this.detector = OverlapDetector.create(intervals);
        this.mode = mode;
    }

    @Override
    public void writeHeader(final VCFHeader header) {
        underlyingWriter.writeHeader(header);
    }

    @Override
    public void setHeader(final VCFHeader header) {
        underlyingWriter.setHeader(header);
    }

    @Override
    public void close() {
        underlyingWriter.close();
        if (filteredCount > 0) {
            logger.info("Removed " + filteredCount + " variants from the output according to '"+mode.getName()+"' variant interval filtering rule.");
        }
    }

    @Override
    public boolean checkError() {
        return underlyingWriter.checkError();
    }

    /**
     * Add the given variant to the writer and output it if it matches.
     * @param vc the variant to potentially write
     */
    @Override
    public void add(final VariantContext vc) {
        if(mode.test(detector, vc)) {
            underlyingWriter.add(vc);
        } else {
            filteredCount++;
        }
    }

}
