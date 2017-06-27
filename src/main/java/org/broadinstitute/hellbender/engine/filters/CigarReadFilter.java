package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collection;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import java.util.ArrayList;

/**
 * Read filter that keeps reads based on a given descriptor for a cigar string.
 * This class utilizes its own syntax for describing a cigar string.
 *
 * Created by jonn on 6/8/17.
 */
public class CigarReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1L;

    //\*|([0-9]+[H])?([0-9]+[S])?([0-9]+[MIDNPX=])+([0-9]+[S])?([0-9]+[H])?
    private static final Pattern validFilterPattern = Pattern.compile (
            "^\\*$|^\\^?(?:\\d*H)*(?:\\d*S)*(?:\\d*[MIDNPX=%])*(?:\\d*S)*(?:\\d*H)*\\$?$"
    );

    private static final Pattern nextCigarMatchElementPattern = Pattern.compile (
            "(\\^?\\d*[MIDNPX=%]\\$?)"
    );

    private String description;

    // =======================================================

    /**
     * Class to hold a single Cigar character (with modifiers).
     */
    private static class CigarMatchElement {

        private CigarOperator operator          = null;

        private int length                      = -1;

        private boolean lessThan                = false;
        private boolean lessThanEqualTo         = false;
        private boolean greaterThan             = false;
        private boolean greaterThanEqualTo      = false;

        private boolean anchoredStart           = false;
        private boolean anchoredEnd             = false;

        private boolean isWildCard              = false;
        private boolean isUnavailable           = false;

        public CigarMatchElement() {}

        public CigarMatchElement(final CigarOperator operator)
        {
            this.operator = operator;
        }

        public CigarOperator getOperator() {
            return operator;
        }

        public void setOperator(CigarOperator operator) {
            this.operator = operator;
        }

        public int getLength() {
            return length;
        }

        public void setLength(int length) {
            this.length = length;
        }

        public boolean isLessThan() {
            return lessThan;
        }

        public void setLessThan(boolean lessThan) {
            this.lessThan = lessThan;
        }

        public boolean isLessThanEqualTo() {
            return lessThanEqualTo;
        }

        public void setLessThanEqualTo(boolean lessThanEqualTo) {
            this.lessThanEqualTo = lessThanEqualTo;
        }

        public boolean isGreaterThan() {
            return greaterThan;
        }

        public void setGreaterThan(boolean greaterThan) {
            this.greaterThan = greaterThan;
        }

        public boolean isGreaterThanEqualTo() {
            return greaterThanEqualTo;
        }

        public void setGreaterThanEqualTo(boolean greaterThanEqualTo) {
            this.greaterThanEqualTo = greaterThanEqualTo;
        }

        public boolean isAnchoredStart() {
            return anchoredStart;
        }

        public void setAnchoredStart(boolean anchoredStart) {
            this.anchoredStart = anchoredStart;
        }

        public boolean isAnchoredEnd() {
            return anchoredEnd;
        }

        public void setAnchoredEnd(boolean anchoredEnd) {
            this.anchoredEnd = anchoredEnd;
        }

        public boolean isWildCard() {
            return isWildCard;
        }

        public void setWildCard(boolean wildCard) {
            isWildCard = wildCard;
        }

        public boolean isUnavailable() {
            return isUnavailable;
        }

        public void setUnavailable(boolean unavailable) {
            isUnavailable = unavailable;
        }
    }

    private Collection<CigarMatchElement> matchElementCollection;

    // =======================================================

    /**
     * Create a {@link CigarReadFilter}.
     */
    public CigarReadFilter() {
        this.setDescription(description);
    }

    /**
     * Create a {@link CigarReadFilter} with the given description.
     * @param description String describing the format of the cigar string on which to filter.
     */
    public CigarReadFilter(final String description) {
        this.setDescription(description);
    }

    public String getDescription() { return description; }

    public void setDescription(final String description) {

        // Before we set our description string, we must make sure it's valid
        if (!validate(description))
        {
            throw new UserException.BadInput("The given cigar filter string is not a valid filter: " + description);
        }

        // Set the description string:
        this.description = description;

        // Parse the string and extract the elements to our list:
        this.matchElementCollection = extractMatchElementsFromString(description);
    }

    /**
     * Parses the description string into individual {@link CigarMatchElement} objects.
     * NOTE: Assumes that the passed description string is valid.
     * @param description String describing the format of the cigar string on which to filter.
     * @return A list of {@link CigarMatchElement} that corresponds to the given description.
     */
    private Collection<CigarMatchElement> extractMatchElementsFromString(final String description) {

        ArrayList<CigarMatchElement> matchElementList = new ArrayList<>();

        // Check for trivial cases here:
        if ( description.compareTo("*") == 0 )
        {
            final CigarMatchElement e = new CigarMatchElement();
            e.setUnavailable(true);
            matchElementList.add(e);
        }
        else
        {
            // Iterate through the string and parse out each match element
            final Matcher matcher = nextCigarMatchElementPattern.matcher(description);
            while (matcher.find())
            {
                // Add a new element to our list based on the string match
                matchElementList.add(createMatchElementFromMatchString(matcher.group()));
            }
        }

        return matchElementList;
    }

    /**
     * Creates a CigarMatchElement from the given string.
     * Assumes the string is a valid {@link CigarMatchElement}.
     * @param match The string representation of a {@link CigarMatchElement}.
     * @return The {@link CigarMatchElement} representation of the given string.
     */
    private CigarMatchElement createMatchElementFromMatchString(final String match) {
        CigarMatchElement e = new CigarMatchElement();

        //TODO: Finish ME!

        // We scan through the 

        return e;
    }

    /**
     * Check whether {@link CigarMatchElement#description} is a valid filter.
     * @return True if {@link CigarMatchElement#description} is valid; False otherwise.
     */
    public boolean validate() {
        return validate(description);
    }

    /**
     * Check whether the given description string is a valid filter.
     * @param description The filter description to validate.
     * @return True if {@code description} is valid; False otherwise.
     */
    public static boolean validate(final String description) {

//     Cigar strings take the (regex) form:
//
//          \^*$|([0-9]+[H])?([0-9]+[S])?([0-9]+[MIDNPX=])+([0-9]+[S])?([0-9]+[H])?
//
//        The significance of each character is the following:
//
//          M alignment match (can be a sequence match or mismatch)
//          I insertion to the reference
//          D deletion from the reference
//          N skipped region from the reference
//          S soft clipping (clipped sequences present in SEQ)
//          H hard clipping (clipped sequences NOT present in SEQ)
//          P padding (silent deletion from padded reference)
//          = sequence match
//          X sequence mismatch
//
//     (See https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf for more details)

        // =======================================================

        // We check that the length is not zero so that we can ensure that if we are on the right-hand side of the match
        // that there is some content there.
        boolean isValid = description.length() != 0;

        if (description.length() == 1) {
            // Since we're allowing for the user to input ^ and $, we need to make sure
            // that the pattern we match against
            isValid = isValid && (description.compareTo("^") != 0);
            isValid = isValid && (description.compareTo("$") != 0);
        }

        isValid = isValid && validFilterPattern.matcher(description).matches();

        return isValid;
    }

    @Override
    public boolean test(GATKRead read) {

        //TODO: Go through each cigar read element and see if it matches.

        return false;
    }
}
