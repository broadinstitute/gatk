package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.List;
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

    /** Regular Expression {@link Pattern} to validate the {@link CigarReadFilter} pattern as given to this {@link CigarReadFilter}. */
    static final Pattern validFilterPattern = Pattern.compile (
            "^\\*$" +
            "|" +
            "^\\^?" +
            "(?:(?:[<>]\\d+|[<>]=\\d+)?\\d*H)*" +
            "(?:(?:[<>]\\d+|[<>]=\\d+)?\\d*S)*" +
            "(?:(?:[<>]\\d+|[<>]=\\d+)?\\d*[MIDNPX=" + CigarMatchElement.WILDCARD + "])*" +
            "(?:(?:[<>]\\d+|[<>]=\\d+)?\\d*S)*" +
            "(?:(?:[<>]\\d+|[<>]=\\d+)?\\d*H)*" +
            "\\$?$"
    );

    /** Regular Expression {@link Pattern} to match the next {@link CigarMatchElement} in the description string. */
    static final Pattern nextCigarMatchElementPattern = Pattern.compile (
            "(\\^?(?:[<>]\\d+|[<>]=\\d+|\\d+)?[SHMIDNPX=" + CigarMatchElement.WILDCARD + "]\\$?)"
    );

    @Argument(fullName = "description",
            shortName = "d",
            doc = "Description corresponding to the kind of cigar strings to include in the filter.",
            optional = true,
            maxElements = 1)
    // Set the default description to be completely permissive.
    private String description = Character.toString( CigarMatchElement.WILDCARD );

    @Argument(fullName = "greedyWildcard",
            doc = "",
            optional = true,
            maxElements = 1)
    private boolean greedyWildcard = false;

    /**
     * The list of CigarMatchElements comprising this {@link CigarReadFilter}.
     */
    private List<CigarMatchElement> matchElementList;

    /**
     * Whether ${@link CigarReadFilter#matchElementList} has an element containing an end anchor
     */
    private boolean hasEndAnchor = false;

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

    public boolean isGreedyWildcard() {
        return greedyWildcard;
    }

    public void setGreedyWildcard(boolean greedyWildcard) {
        this.greedyWildcard = greedyWildcard;
    }

    public String getDescription() { return description; }

    /**
     * Validates the given description and then stores it as {@link CigarReadFilter#description}.
     * @param description String describing the format of the cigar string on which to filter.
     */
    public void setDescription(final String description) {

        // Before we set our description string, we must make sure it's valid
        if (!validate(description))
        {
            throw new UserException.BadInput("The given cigar filter string is not a valid filter: " + description);
        }

        // Set the description string:
        this.description = description;

        // Set our end anchor information:
        this.hasEndAnchor = description.contains("$");

        // Parse the string and extract the elements to our list:
        this.matchElementList = extractMatchElementsFromString(description);
    }

    /**
     * Parses the description string into individual {@link CigarMatchElement} objects.
     * NOTE: Assumes that the passed description string is valid.
     * @param description String describing the format of the cigar string on which to filter.
     * @return A list of {@link CigarMatchElement} that corresponds to the given description.
     */
    List<CigarMatchElement> extractMatchElementsFromString(final String description) {

        ArrayList<CigarMatchElement> matchElementList = new ArrayList<>();

        // Check for trivial cases here:
        if ( description.compareTo(CigarMatchElement.UNAVAILABLE) == 0 )
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
                matchElementList.add(createMatchElementFromSingleMatchString(matcher.group()));
            }
        }

        return matchElementList;
    }

    /**
     * Creates a CigarMatchElement from the given string.
     * Assumes the string is a valid representation of a {@link CigarMatchElement}
     * and that it does not code for an unavailable Cigar String.
     * @param match The string representation of a {@link CigarMatchElement}.
     * @return The {@link CigarMatchElement} representation of the given string.
     */
    CigarMatchElement createMatchElementFromSingleMatchString(String match) {
        CigarMatchElement e = new CigarMatchElement();

        // If we have a ^ character, we set the appropriate flag and cut it off the front:
        if ( match.charAt(0) == '^' ) {
            e.setAnchoredStart(true);
            match = match.substring(1);
        }

        // If we have a $ character, we set the appropriate flag and cut it off the end:
        if ( match.charAt(match.length()-1) == '$' ) {
            e.setAnchoredEnd(true);
            match = match.substring(0, match.length()-1);
        }

        // If we have qualifiers for numerical values:
        // NOTE: We know these must be at the new start of the string if they exist,
        //       so we can remove them from the start of the string if they're there.
        if ( match.charAt(0) == '<' ) {
            if ( match.charAt(1) == '=' ) {
                e.setLessThanEqualTo(true);
                match = match.substring(2);
            }
            else {
                e.setLessThan(true);
                match = match.substring(1);
            }
        }
        else if ( match.charAt(0) == '>' ) {
            if ( match.charAt(1) == '=' ) {
                e.setGreaterThanEqualTo(true);
                match = match.substring(2);
            }
            else {
                e.setGreaterThan(true);
                match = match.substring(1);
            }
        }

        // At this point, all that's left is a numerical designator at the front of the string
        // and a character representing the cigar type.

        // Get the cigar type now:
        if ( match.charAt(match.length()-1) == CigarMatchElement.WILDCARD ) {
            e.setWildCard(true);
        }
        else {
            e.setOperator( CigarOperator.characterToEnum( match.charAt(match.length()-1)) );
        }
        match = match.substring(0, match.length()-1);

        // Now we can just use the rest of the string as a number:
        if (match.length() > 0) {
            e.setLength( Integer.valueOf(match) );
        }

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
//        \\*$|^(?:\\d*H)?(?:\\d*S)?(?:\\d*[MIDNPX=%])*(?:\\d*S)*(?:\\d*H)*$
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
            // that the pattern we match against doesn't consist of only those anchors
            isValid = isValid && (description.compareTo("^") != 0);
            isValid = isValid && (description.compareTo("$") != 0);
        }

        isValid = isValid && validFilterPattern.matcher(description).matches();

        return isValid;
    }

    /**
     * Tests the given {@link GATKRead} to see if its cigar string matches the specified match description in {@link CigarMatchElement#description}.
     * NOTE: Assumes the given read's cigar string and the match description {@link CigarMatchElement#description} are both valid.
     * @param read The {@link GATKRead} containing the {@link Cigar} to check for inclusion in this {@link CigarReadFilter}
     * @return {@code true} if the given read's cigar string matches {@link CigarMatchElement#description}, {@code false} otherwise.
     */
    @Override
    public boolean test(GATKRead read) {

        //We'll need to check each cigar element and can't assume a normalized cigar.
        //This is OK.  It'll be faster to do this than to normalize and then check it.

        Cigar cigar = read.getCigar();

        // Quick check of the cigar to see if it's unavailable and if we're looking for the unavailable case:
        if (cigar.toString() == "*") {
            if ((matchElementList.size() == 1) && (matchElementList.get(0).isUnavailable())) {
                return true;
            }
            else {
                return false;
            }
        }

        // The tests are a little different if we are using greedy wildcards or not.
        // Make sure we use the right one.
        if ( greedyWildcard ) {
            return testWithGreedyWildcard(cigar);

        }
        else {
            return testWithoutGreedyWildcard(cigar);
        }
    }

    /**
     * Tests the given {@link Cigar} to see if its cigar string matches the specified match description in {@link CigarMatchElement#description}.
     * NOTE: Assumes the given read's cigar string and the match description {@link CigarMatchElement#description} are both valid
     *       and that the matching of any wildcards is to be done in a non-greedy fashion.
     * @param cigar The {@link Cigar} to check for inclusion in this {@link CigarReadFilter}
     * @return {@code true} if the given read's cigar string matches {@link CigarMatchElement#description}, {@code false} otherwise.
     */
    private boolean testWithoutGreedyWildcard(Cigar cigar) {

        // Keep track of our position in the cigar string out here, so we can have control over where we are
        int cigarIndex = 0;

        // If we have an anchored match element, then we need to start from
        // the point in the cigar string at which the match can begin:
        if ( hasEndAnchor ) {
            // The match element of the end anchor must be the last element.
            // Therefore we do some simple math to get to the start of where the match
            // should begin:
            cigarIndex = cigar.numCigarElements() - matchElementList.size();

            // check for degenerate case:
            if (cigarIndex < 0) {
                return false;
            }
        }

        // Keep track of if we've seen more than one operator (for anchor checking):
        boolean haveSeenMoreThanOneOperator = false;

        // Keep track of whether our match has begun:
        boolean hasBegunMatch = false;

        // We go through our match elements to make sure that the cigar string matches each one:
        for ( CigarMatchElement matchElement : matchElementList ) {

            // Make sure we can keep track of the number of operators we've seen:
            int cigarOperatorCount = 0;
            CigarElement lastElementProcessed = null;

            // Check if we have more cigar string to check.
            // If we do not, then we don't match:
            if ( cigarIndex == cigar.numCigarElements() ) {
                // We have no more cigar elements to check, but we still
                // have a match element to check against.
                // This indicates failure.
                return false;
            }

            for ( ; cigarIndex < cigar.numCigarElements() ; ++cigarIndex ) {

                // Get the next cigar element:
                CigarElement cigarElement = cigar.getCigarElement(cigarIndex);

                // We have a new element.  Do some checks involving when the operator changes:
                if ((lastElementProcessed != null) && (lastElementProcessed.getOperator().ordinal() != cigarElement.getOperator().ordinal())) {

                    // Update our flag for more than one operator:
                    haveSeenMoreThanOneOperator = true;

                    // If we have gotten here, then we're still OK
                    // but we need to get the next matchElement, so we break
                    // to cause the outer loop to iterate:
                    break;
                }

                // Check to see if the operator must be anchored:
                if ( matchElement.isAnchoredStart() ) {
                    // If we're an anchored start, we only want to continue
                    // if we have only seen one operator since the start of the string:
                    if ( haveSeenMoreThanOneOperator ) {
                        return false;
                    }
                }

                // Check the cigar operator first:
                if ( (!matchElement.isWildCard()) &&
                        (cigarElement.getOperator().ordinal() != matchElement.getOperator().ordinal()) ) {

                    // Our cigar element and our match element aren't the same.

                    // If we have not started matching yet, we advance to the next cigar element:
                    if (!hasBegunMatch) {
                        continue;
                    }
                    else {
                        // We have already begun matching and our element doesn't match.
                        // This is a failure case.
                        return false;
                    }
                }

                // Increment the length of the operators so far:
                cigarOperatorCount += cigarElement.getLength();

                // If we need to worry about length:
                if ( matchElement.requiresLength() ) {

                    // Check the "less than" cases here at the end of the loop:
                    if (matchElement.isLessThan() && (cigarOperatorCount >= matchElement.getLength())) {
                        // We have too many elements of this type.
                        // We don't match.
                        return false;
                    } else if (matchElement.isLessThanEqualTo() && (cigarOperatorCount > matchElement.getLength())) {
                        // We have too many elements of this type.
                        // We don't match.
                        return false;
                    }
                }

                // Store the last element that we looked at so we can do some comparison work later:
                lastElementProcessed = cigarElement;

                // We have started matching:
                hasBegunMatch = true;
            }

            // Now check the "greater-than" and "equals" length cases if this is a new (different) element than the last one
            // and only if we care about length:
            if (matchElement.requiresLength()) {
                if (matchElement.isGreaterThan()) {
                    if (cigarOperatorCount <= matchElement.getLength()) {
                        // We have too few operators of the right type.
                        // We do not pass:
                        return false;
                    }
                } else if (matchElement.isGreaterThanEqualTo() ) {
                    if (cigarOperatorCount < matchElement.getLength()) {
                        // We have too few operators of the right type.
                        // We do not pass:
                        return false;
                    }
                    //The strictly equals case:
                } else if ( !(matchElement.isLessThanEqualTo() || matchElement.isLessThan()) ) {
                    if (cigarOperatorCount != matchElement.getLength()) {
                        // We do not have the exact number of operators required.
                        // We do not pass:
                        return false;
                    }
                }
            }
        }

        // If we get out of all the checks and have gotten to the end
        // of our cigar string and our match string, then if we have matched
        // at all, the whole string is a match:
        return hasBegunMatch;
    }

    /**
     * Tests the given {@link Cigar} to see if its cigar string matches the specified match description in {@link CigarMatchElement#description}.
     * NOTE: Assumes the given read's cigar string and the match description {@link CigarMatchElement#description} are both valid
     *       and that the matching of any wildcards is to be done in a greedy fashion.
     * @param rawCigar The {@link Cigar} to check for inclusion in this {@link CigarReadFilter}
     * @return {@code true} if the given read's cigar string matches {@link CigarMatchElement#description}, {@code false} otherwise.
     */
    private boolean testWithGreedyWildcard(Cigar rawCigar) {
        // TODO: FIXME

        return false;

//        Cigar cigar = CigarUtils.combineAdjacentCigarElements( rawCigar );
//
//        // Keep track of our position in the cigar string out here, so we can have control over where we are
//        int cigarIndex = 0;
//
//        // If we have an anchored match element, then we need to start from
//        // the point in the cigar string at which the match can begin:
//        if ( hasEndAnchor ) {
//            // The match element of the end anchor must be the last element.
//            // Therefore we do some simple math to get to the start of where the match
//            // should begin:
//            cigarIndex = cigar.numCigarElements() - matchElementList.size();
//
//            // check for degenerate case:
//            if (cigarIndex < 0) {
//                return false;
//            }
//        }
//
//        // Keep track of if we've seen more than one operator (for anchor checking):
//        boolean haveSeenMoreThanOneOperator = false;
//
//        // Keep track of whether our match has begun:
//        boolean hasBegunMatch = false;
//
//        // A cigar element to keep track of remaining parts of a previous element.
//        // We need this for the wildcard matches.
//        CigarElement remainderElement = null;
//
//        // We go through our match elements to make sure that the cigar string matches each one:
//        for ( int matchElementIndex = 0; matchElementIndex < matchElementList.size() ; ++matchElementIndex ) {
//
//            CigarMatchElement matchElement = matchElementList.get(matchElementIndex);
//
//            // Peek at the next cigar match element:
//            CigarMatchElement nextMatch = (matchElementIndex < matchElementList.size() - 1)
//                    ? matchElementList.get(matchElementIndex + 1)
//                    : null;
//
//            // Make sure we can keep track of the number of operators we've seen:
//            int cigarOperatorCount = 0;
//            int wildcardOperatorCount = 0;
//            CigarElement lastElementProcessed = null;
//
//            // Check if we have more cigar string to check.
//            // If we do not, then we don't match:
//            if ( cigarIndex == cigar.numCigarElements() ) {
//                // We have no more cigar elements to check, but we still
//                // have a match element to check against.
//                // This indicates failure.
//                return false;
//            }
//
//            while ( cigarIndex < cigar.numCigarElements() ) {
//
//                if ( remainderElement == null ) {
//                    // Get the next cigar element:
//                    remainderElement = cigar.getCigarElement(cigarIndex);
//                    ++cigarIndex;
//                }
//
//                // Special handling for wildcards:
//                if ( matchElement.isWildCard() ) {
//
//                    // Do we have a numerical qualifier?
//                    if ( matchElement.isLessThanEqualTo() ) {
//                        int wildcardOperatorsNeeded = matchElement.getLength() - wildcardOperatorCount;
//
//                        if ( remainderElement.getLength() < wildcardOperatorsNeeded ) {
//                            // We don't have enough elements to fill this wildcard.
//                            // Add these to the wildcard count and continue to the next CigarElement
//                            wildcardOperatorCount += matchElement.getLength();
//                            remainderElement = null;
//                            continue;
//                        }
//                        else if ( remainderElement.getLength() > wildcardOperatorsNeeded ) {
//                            // This match element has too many elements for us.
//                            // We need to subtract the length and hold the remaining elements in our
//                            // remainder object.
//                            // Then we need to try this loop again with the remainder element:
//                            remainderElement = new CigarElement(
//                                    remainderElement.getLength() - wildcardOperatorsNeeded,
//                                    remainderElement.getOperator()
//                            );
//                            continue;
//                        }
//                        else { // Equals case.
//                            // We have exactly the number of operators needed for this wildcard.
//                            // Go to the next CigarMatchElement and the next CigarElement
//                            remainderElement = null;
//                            break;
//                        }
//                    }
//                    else if ( matchElement.isLessThan() ) {
//
//                    }
//                    else if ( matchElement.isGreaterThanEqualTo() ) {
//
//                    }
//                    else if ( matchElement.isGreaterThan() ) {
//
//                    }
//                    else {
//
//                        if ( nextMatch == null ) {
//                            // We have no more match elements and we are a wildcard.
//                            // This means we trivially match to the end of the string.
//                            return true;
//                        }
//                        else
//                        {
//                            // We have another match element, so we need to determine when we can break out of the
//                            // wildcard matching.
//
//                            // keep going until we find a cigar element that matches the next match element
//                        }
//
//                    }
//                }
//
//                // We have a new element.  Do some checks involving when the operator changes:
//                if ((lastElementProcessed != null) && (lastElementProcessed.getOperator().ordinal() != cigarElement.getOperator().ordinal())) {
//
//                    // Update our flag for more than one operator:
//                    haveSeenMoreThanOneOperator = true;
//
//                    // If we have gotten here, then we're still OK
//                    // but we need to get the next matchElement, so we break
//                    // to cause the outer loop to iterate:
//                    break;
//                }
//
//                // Check to see if the operator must be anchored:
//                if ( matchElement.isAnchoredStart() ) {
//                    // If we're an anchored start, we only want to continue
//                    // if we have only seen one operator since the start of the string:
//                    if ( haveSeenMoreThanOneOperator ) {
//                        return false;
//                    }
//                }
//
//                // Check the cigar operator first:
//                if ( (!matchElement.isWildCard()) &&
//                        (cigarElement.getOperator().ordinal() != matchElement.getOperator().ordinal()) ) {
//
//                    // Our cigar element and our match element aren't the same.
//
//                    // If we have not started matching yet, we advance to the next cigar element:
//                    if (!hasBegunMatch) {
//                        ++cigarIndex;
//                        continue;
//                    }
//                    else {
//                        // We have already begun matching and our element doesn't match.
//                        // This is a failure case.
//                        return false;
//                    }
//                }
//
//                // Increment the length of the operators so far:
//                cigarOperatorCount += cigarElement.getLength();
//
//                // If we need to worry about length:
//                if ( matchElement.requiresLength() ) {
//
//                    // Check the "less than" cases here at the end of the loop:
//                    if (matchElement.isLessThan() && (cigarOperatorCount >= matchElement.getLength())) {
//                        // We have too many elements of this type.
//                        // We don't match.
//                        return false;
//                    } else if (matchElement.isLessThanEqualTo() && (cigarOperatorCount > matchElement.getLength())) {
//                        // We have too many elements of this type.
//                        // We don't match.
//                        return false;
//                    }
//                }
//
//                // Store the last element that we looked at so we can do some comparison work later:
//                lastElementProcessed = cigarElement;
//
//                // We have started matching:
//                hasBegunMatch = true;
//
//                ++cigarIndex;
//            }
//
//            // Now check the "greater-than" and "equals" length cases if this is a new (different) element than the last one
//            // and only if we care about length:
//            if (matchElement.requiresLength()) {
//                if (matchElement.isGreaterThan()) {
//                    if (cigarOperatorCount <= matchElement.getLength()) {
//                        // We have too few operators of the right type.
//                        // We do not pass:
//                        return false;
//                    }
//                } else if (matchElement.isGreaterThanEqualTo() ) {
//                    if (cigarOperatorCount < matchElement.getLength()) {
//                        // We have too few operators of the right type.
//                        // We do not pass:
//                        return false;
//                    }
//                    //The strictly equals case:
//                } else if ( !(matchElement.isLessThanEqualTo() || matchElement.isLessThan()) ) {
//                    if (cigarOperatorCount != matchElement.getLength()) {
//                        // We do not have the exact number of operators required.
//                        // We do not pass:
//                        return false;
//                    }
//                }
//            }
//        }
//
//        // If we get out of all the checks and have gotten to the end
//        // of our cigar string and our match string, then if we have matched
//        // at all, the whole string is a match:
//        return hasBegunMatch;
    }

    // =======================================================

    /**
     * Class to hold a single Cigar character for {@link CigarReadFilter#description} (with modifiers).
     * This is different from a {@link CigarElement}, which holds a single
     * operator for a Cigar string.
     */
    static class CigarMatchElement {

        public static final String UNAVAILABLE  = "*";
        public static final char WILDCARD       = '%';

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

        public boolean requiresLength() { return length != -1; }

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

        @Override
        public boolean equals(Object o) {
            if ( o == this ) {
                return true;
            }

            if ( !(o instanceof CigarMatchElement) ) {
                return false;
            }

            CigarMatchElement c = (CigarMatchElement)o;

            return  (c.operator           == this.operator)           &&
                    (c.length             == this.length)             &&
                    (c.lessThan           == this.lessThan)           &&
                    (c.lessThanEqualTo    == this.lessThanEqualTo)    &&
                    (c.greaterThan        == this.greaterThan)        &&
                    (c.greaterThanEqualTo == this.greaterThanEqualTo) &&
                    (c.anchoredStart      == this.anchoredStart)      &&
                    (c.anchoredEnd        == this.anchoredEnd)        &&
                    (c.isWildCard         == this.isWildCard)         &&
                    (c.isUnavailable      == this.isUnavailable);
        }

        @Override
        public int hashCode() {
            StringBuilder stringBuilder = new StringBuilder();

            stringBuilder.append(this.operator);
            stringBuilder.append(this.length);
            stringBuilder.append(this.lessThan);
            stringBuilder.append(this.lessThanEqualTo);
            stringBuilder.append(this.greaterThan);
            stringBuilder.append(this.greaterThanEqualTo);
            stringBuilder.append(this.anchoredStart);
            stringBuilder.append(this.anchoredEnd);
            stringBuilder.append(this.isWildCard);
            stringBuilder.append(this.isUnavailable);

            return stringBuilder.toString().hashCode();
        }
    }
}
