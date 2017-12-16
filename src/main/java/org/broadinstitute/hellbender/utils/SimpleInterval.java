 package org.broadinstitute.hellbender.utils;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.mergealignment.SamAlignmentMerger;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

 /**
 * Minimal immutable class representing a 1-based closed ended genomic interval
 * SimpleInterval does not allow null contig names.  It cannot represent an unmapped Locatable.
 *
 *@warning 0 length intervals are NOT currently allowed, but support may be added in the future
 */
public final class SimpleInterval implements Locatable, Serializable {
     private static final Logger log = LogManager.getLogger(SimpleInterval.class);

    private static final long serialVersionUID = 1L;
    public static final char CONTIG_SEPARATOR = ':';
    public static final char START_END_SEPARATOR = '-';
    public static final String END_OF_CONTIG = "+"; //note: needs to be a String because it's used in an endsWith call.

    private final int start;
    private final int end;
    private final String contig;

    /**
     * Create a new immutable 1-based interval of the form [start, end]
     * @param contig the name of the contig, must not be null
     * @param start  1-based inclusive start position
     * @param end  1-based inclusive end position
     */
    public SimpleInterval(final String contig, final int start, final int end){
        validatePositions(contig, start, end);
        this.contig = contig;
        this.start = start;
        this.end = end;
    }

    /**
     * Create a new SimpleInterval from a {@link Locatable}
     * @param locatable any Locatable
     * @throws IllegalArgumentException if locatable violates any of the SimpleInterval constraints or is null
     */
    public SimpleInterval(final Locatable locatable){
        this(Utils.nonNull(locatable).getContig(),
                locatable.getStart(), locatable.getEnd());
    }

    /**
     * Test that these are valid values for constructing a SimpleInterval:
     *    contig cannot be null
     *    start must be >= 1
     *    end must be >= start
     * @throws IllegalArgumentException if it is invalid
     */
    private static void validatePositions(final String contig, final int start, final int end) {
        Utils.validateArg(isValid(contig, start, end), () -> "Invalid interval. Contig:" + contig + " start:"+start + " end:" + end);
    }

     /**
      * Test that these are valid values for constructing a SimpleInterval:
      *    contig cannot be null
      *    start must be >= 1
      *    end must be >= start
      */
     public static boolean isValid(final String contig, final int start, final int end) {
         return contig != null && start > 0 && end >= start;
     }

     /**
     * Makes an interval by parsing the string.
     *
     * @warning this method does not fill in the true contig end values
     * for intervals that reach to the end of their contig,
     * uses {@link Integer#MAX_VALUE} instead.
     *
     * Semantics of start and end are defined in {@link Locatable}.
     * The format is one of:
     *
     * contig           (Represents the whole contig, from position 1 to the {@link Integer#MAX_VALUE})
     *
     * contig:start     (Represents the 1-element range start-start on the given contig)
     *
     * contig:start-end (Represents the range start-end on the given contig)
     *
     * contig:start+    (Represents the prefix of the contig starting at the given start position and ending at {@link Integer#MAX_VALUE})
     *
     * examples (note that _all_ commas in numbers are simply ignored, for human convenience):
     *
     * 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000' or 'chr2:1000000+'
      *
      * @param str non-empty string to be parsed
     */
    public SimpleInterval(final String str){
        /* Note: we want to keep the class immutable. So all fields need to be final.
         * But only constructors can assign to final fields.
         * So we can either keep this parsing code in the constructor or make a static factory method
         * and make multiple objects. We chose the former.
         */
        Utils.nonNull(str);
        Utils.validateArg(!str.isEmpty(), "str should not be empty");

        final String contig;
        final int start;
        final int end;

        final int colonIndex = str.lastIndexOf(CONTIG_SEPARATOR);
        if (colonIndex == -1) {
            contig = str;  // chr1
            start = 1;
            end = Integer.MAX_VALUE;
        } else {
            contig = str.substring(0, colonIndex);
            final int dashIndex = str.indexOf(START_END_SEPARATOR, colonIndex);
            if(dashIndex == -1) {
                if(str.endsWith(END_OF_CONTIG)) {
                    start = parsePosition(str.substring(colonIndex + 1, str.length() - 1));  // chr:1+
                    end = Integer.MAX_VALUE;
                } else {
                    start = parsePosition(str.substring(colonIndex + 1));   // chr1:1
                    end = start;
                }
            } else {
                start = parsePosition(str.substring(colonIndex + 1, dashIndex));  // chr1:1-1
                end = parsePosition(str.substring(dashIndex + 1));
            }
        }

        validatePositions(contig, start, end);
        this.contig = contig;
        this.start = start;
        this.end = end;
    }

     /**
      * Given an interval query string and a sequence dictionary, determine if the query string can be
      * resolved as a valid interval query against more than one contig in the dictionary, i.e., more than
      * one of:
      *
      *     prefix
      *     prefix:nnn
      *     prefix:nnn+
      *     prefix:nnn-nnn
      *
      * and return the list of all possible interpretations (there can never be more than 2). Note that for
      * an ambiguity to exist, the query string must contain at least one colon.
      *
      * @param intervalQueryString
      * @param sequenceDictionary
      * @return List<SimpleInterval> containing 0, 1, or 2 valid interpretations of {code queryString} given
      * {@code sequenceDictionary}. If the list is empty, the query doesn't match any contig in the sequence
      * dictionary. If the list contains more than one interval, the query string is ambiguous and should be
      * rejected. If the list contains a single interval, the query is unambiguous and can be safely used to
      * conduct a query.
      * @thows IllegalArgumentException if the query only matches a single contig in the dictionary, but the start
      * and end positions are not valis
      * @throws NumberFormatException if the query only matches a single contig in the dictionary, but the query
      * interval paramaters (start, end) cannot be parsed
      */
     public static List<SimpleInterval> getResolvedIntervals(
             final String intervalQueryString,
             final SAMSequenceDictionary sequenceDictionary) {
         Utils.nonNull(intervalQueryString);
         Utils.validateArg(!intervalQueryString.isEmpty(), "intervalQueryString should not be empty");

         // Keep a list of all valid interpretations
         final List<SimpleInterval> resolvedIntervals = new ArrayList<>();

         // Treat the entire query string as a contig name. If it exists in the sequence dictionary,
         // count that as one valid interpretation.
         final SAMSequenceRecord queryAsContigName = sequenceDictionary.getSequence(intervalQueryString);
         if (queryAsContigName != null) {
             resolvedIntervals.add(new SimpleInterval(intervalQueryString, 1, queryAsContigName.getSequenceLength()));
         }

         // The query must contain at least one colon for an ambiguity to exist.
         final int lastColonIndex = intervalQueryString.lastIndexOf(CONTIG_SEPARATOR);
         if (lastColonIndex == -1) {
             return resolvedIntervals;
         }

         // Get a prefix containing everything up to the last colon, and see if it represents a valid contig.
         final String prefix = intervalQueryString.substring(0, lastColonIndex);
         final SAMSequenceRecord prefixSequence = sequenceDictionary.getSequence(prefix);
         if (prefixSequence == null) {
             return resolvedIntervals;
         }

         try {
             final int lastDashIndex = intervalQueryString.lastIndexOf(START_END_SEPARATOR);
             int startPos;
             int endPos;

             // Try to resolve the suffix as a query against the contig represented by the prefix.
             if (intervalQueryString.endsWith(END_OF_CONTIG)) {
                 // try to resolve as "prefix:nnn+"
                 startPos = parsePositionThrowOnFailure(intervalQueryString.substring(lastColonIndex + 1, intervalQueryString.length()-1));
                 endPos = prefixSequence.getSequenceLength();
             } else if (lastDashIndex > lastColonIndex) {
                 // try to resolve as "prefix:start-end"
                 startPos = parsePositionThrowOnFailure(intervalQueryString.substring(lastColonIndex + 1, lastDashIndex));
                 endPos = parsePositionThrowOnFailure(intervalQueryString.substring(lastDashIndex + 1, intervalQueryString.length()));
             } else {
                 // finally, try to resolve as "prefix:nnn"
                 startPos = parsePositionThrowOnFailure(intervalQueryString.substring(lastColonIndex + 1, intervalQueryString.length()));
                 endPos = startPos;
             }

             if (isValid(prefix, startPos, endPos)) {
                 // We've pre-tested to validate the positions, so add this interval. This should never throw.
                 resolvedIntervals.add(new SimpleInterval(prefix, startPos, endPos));
             } else {
                 // Positions don't appear to be valid, but we don't want to throw if there is any other valid
                 // interpretation of the query string
                 if (resolvedIntervals.isEmpty()) {
                     // validatePositions throws on validation failure, which is guaranteed if we got here
                     validatePositions(prefix, startPos, endPos);
                 }
             }
         } catch (NumberFormatException e) {
             // parsing of the start or end pos failed
             if (resolvedIntervals.isEmpty()) {
                 throw e;
             } else {
                 // We're interpreting this as a query against a full contig, but its POSSIBLE that the user
                 // mis-entered the start or stop position, and had they entered them correctly, would have resulted
                 // in an ambiguity. So accept the query as an interval for the full contig, but issue a warning
                 // saying how the query as resolved.
                 log.warn(String.format(
                         "The query interval string \"%s\" is interpreted as a query against the contig named \"%s\", " +
                                 "but may have been intended as an (accidentally malformed) query against the contig named \"%s\"",
                                 intervalQueryString,
                                 resolvedIntervals.get(0).getContig(),
                                 prefixSequence.getSequenceName()));
             }
         }
         return resolvedIntervals;
     }

     /**
      * Parses a number like 100000 or 1,000,000 into an int. Throws NumberFormatException on parse failure.
      */
     private static int parsePositionThrowOnFailure(final String pos) throws NumberFormatException {
         return Integer.parseInt(pos.replaceAll(",", "")); //strip commas
     }

     /**
     * Parses a number like 100000 or 1,000,000 into an int.
     */
    private static int parsePosition(final String pos) {
        try {
            return parsePositionThrowOnFailure(pos);
        } catch (NumberFormatException e){
            throw new UserException("Problem parsing start/end value in interval string. Value was: " + pos, e);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SimpleInterval that = (SimpleInterval) o;

        if (end != that.end) return false;
        if (start != that.start) return false;
        return contig.equals(that.contig);
    }

    @Override
    public int hashCode() {
        int result = start;
        result = 31 * result + end;
        result = 31 * result + contig.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return IntervalUtils.locatableToString(this);
    }

     /**
     * @return name of the contig this is mapped to
     */
    @Override
    public String getContig(){
        return contig;
    }

    /** Gets the 1-based start position of the interval on the contig. */
    @Override
    public int getStart(){
        return start;
    }

    /**
    * @return the 0-based start position (from the GA4GH spec).
    */
    public long getGA4GHStart() {return start - 1; }

    /**
     * @return the 1-based closed-ended end position of the interval on the contig.
     */
    @Override
    public int getEnd(){
        return end;
    }

    /**
    * @return the typical end spans are [zero-start,end) (from the GA4GH spec).
    */
    public long getGA4GHEnd() { return end; }

    /**
     * @return number of bases covered by this interval (will always be > 0)
     */
    public int size() {
        return end - start + 1;
    }

    /**
     * Determines whether this interval overlaps the provided locatable.
     *
     * @param other interval to check
     * @return true if this interval overlaps other, otherwise false
     */
    public boolean overlaps( final Locatable other ) {
        return overlapsWithMargin(other, 0);
    }

     /**
      * Determines whether this interval comes within "margin" of overlapping the provided locatable.
      * This is the same as plain overlaps if margin=0.
      *
      * @param other interval to check
      * @param margin how many bases may be between the two interval for us to still consider them overlapping; must be non-negative
      * @return true if this interval overlaps other, otherwise false
      * @throws IllegalArgumentException if margin is negative
      */
     public boolean overlapsWithMargin(final Locatable other, final int margin) {
         if ( margin < 0 ) {
             throw new IllegalArgumentException("given margin is negative: " + margin +
                     "\tfor this: " + toString() + "\tand that: " + (other == null ? "other is null" : other.toString()));
         }
         if ( other == null || other.getContig() == null ) {
             return false;
         }

         return this.contig.equals(other.getContig()) && this.start <= other.getEnd() + margin && other.getStart() - margin <= this.end;
     }


     /**
     * Determines whether this interval contains the entire region represented by other
     * (in other words, whether it covers it).
     *
     * @param other interval to check
     * @return true if this interval contains all of the bases spanned by other, otherwise false
     */
    public boolean contains( final Locatable other ) {
        if ( other == null || other.getContig() == null ) {
            return false;
        }

        return this.contig.equals(other.getContig()) && this.start <= other.getStart() && this.end >= other.getEnd();
    }

     /**
      * Returns the intersection of the two intervals. The intervals must overlap or IllegalArgumentException will be thrown.
      */
     public SimpleInterval intersect( final Locatable that ) {
         Utils.validateArg(this.overlaps(that), () ->
                 "SimpleInterval::intersect(): The two intervals need to overlap " + this + " " + that);

         return new SimpleInterval(getContig(),
                 Math.max(getStart(), that.getStart()),
                 Math.min( getEnd(), that.getEnd()) );
     }

     /**
      * Returns a new SimpleInterval that represents the entire span of this and that.  Requires that
      * this and that SimpleInterval are contiguous.
      */
     public SimpleInterval mergeWithContiguous( final Locatable that ) {
         Utils.nonNull(that);
         if (!this.contiguous(that)) {
             throw new GATKException("The two intervals need to be contiguous: " + this + " " + that);
         }

         return new SimpleInterval(getContig(),
                 Math.min( getStart(), that.getStart() ),
                 Math.max( getEnd(), that.getEnd()) );
     }

     /**
      * Returns a new SimpleInterval that represents the region between the endpoints of this and other.
      *
      * Unlike {@link #mergeWithContiguous}, the two intervals do not need to be contiguous
      *
      * @param other the other interval with which to calculate the span
      * @return a new SimpleInterval that represents the region between the endpoints of this and other.
      */
     public SimpleInterval spanWith( final Locatable other ) {
         Utils.nonNull(other);
         Utils.validateArg(this.getContig().equals(other.getContig()), "Cannot get span for intervals on different contigs");
         return new SimpleInterval(contig, Math.min(start, other.getStart()), Math.max(end, other.getEnd()));
     }

     private boolean contiguous(final Locatable that) {
         Utils.nonNull(that);
         return this.getContig().equals(that.getContig()) && this.getStart() <= that.getEnd() + 1 && that.getStart() <= this.getEnd() + 1;
     }

     /**
      * Returns a new SimpleInterval that represents this interval as expanded by the specified amount in both
      * directions, bounded by the contig start/stop if necessary.
      *
      * @param padding amount to expand this interval
      * @param contigLength length of this interval's contig
      * @return a new SimpleInterval that represents this interval as expanded by the specified amount in both
      *         directions, bounded by the contig start/stop if necessary.
      */
     public SimpleInterval expandWithinContig( final int padding, final int contigLength ) {
         Utils.validateArg(padding >= 0, "padding must be >= 0");
         return IntervalUtils.trimIntervalToContig(contig, start - padding, end + padding, contigLength);
     }

     /**
      * Returns a new SimpleInterval that represents this interval as expanded by the specified amount in both
      * directions, bounded by the contig start/stop if necessary.
      *
      * @param padding amount to expand this interval
      * @param sequenceDictionary dictionary to use to determine the length of this interval's contig
      * @return a new SimpleInterval that represents this interval as expanded by the specified amount in both
      *         directions, bounded by the contig start/stop if necessary.
      */
     public SimpleInterval expandWithinContig( final int padding, final SAMSequenceDictionary sequenceDictionary ) {
         Utils.nonNull(sequenceDictionary);
         final SAMSequenceRecord contigRecord = sequenceDictionary.getSequence(contig);
         Utils.nonNull( contigRecord, () -> "Contig " + contig + " not found in provided dictionary");

         return expandWithinContig(padding, contigRecord.getSequenceLength());
     }
 }
