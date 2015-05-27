package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.text.XReadLines;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Parse text representations of interval strings that
 * can appear in GATK-based applications.
 */
public final class IntervalUtils {

    /**
     * Lexicographical (contig) order comparator.
     * <p>
     * Intervals from different contigs order is according their enclosing
     * contigs name ascending lexicographical order.
     * </p>
     * <p>
     * Intervals from the same contigs order is according to their start
     * position ascending numerical order, and, in case of a tie, the stop position's.
     * </p>
     * <p>
     * The {@code null} contig is supported and comes last.
     * </p>
     */
    public final static Comparator<Locatable> LEXICOGRAPHICAL_ORDER_COMPARATOR =
            Comparator.comparing(Locatable::getContig,Comparator.nullsLast(String::compareTo))
                    .thenComparingInt(Locatable::getStart)
                    .thenComparingInt(Locatable::getEnd);

    private static Logger logger = LogManager.getLogger(IntervalUtils.class);

    public static GenomeLocSortedSet loadIntervals(
            final List<String> intervalStrings,
            final IntervalSetRule intervalSetRule,
            final IntervalMergingRule intervalMergingRule,
            final int padding,
            final GenomeLocParser genomeLocParser) {
        List<GenomeLoc> allIntervals = new ArrayList<>();
        for ( String intervalString : intervalStrings) {
            List<GenomeLoc> intervals = parseIntervalArguments(genomeLocParser, intervalString);

            if ( padding > 0 ) {
                intervals = getIntervalsWithFlanks(genomeLocParser, intervals, padding);
            }

            allIntervals = mergeListsBySetOperator(intervals, allIntervals, intervalSetRule);
        }

        return sortAndMergeIntervals(genomeLocParser, allIntervals, intervalMergingRule);
    }


    /**
     * Turns a set of strings describing intervals into a parsed set of intervals.  Valid string elements can be files,
     * intervals in samtools notation (chrA:B-C), or some combination of the above separated by semicolons.  Additionally,
     * 'all' can be supplied to indicate all possible intervals, but 'all' must be exclusive of all other interval
     * specifications.
     *
     * @param parser Genome loc parser.
     * @param argList A list of strings containing interval data.
     * @return an unsorted, unmerged representation of the given intervals.  Null is used to indicate that all intervals should be used.
     */
    public static List<GenomeLoc> parseIntervalArguments(GenomeLocParser parser, List<String> argList) {
        List<GenomeLoc> rawIntervals = new ArrayList<>();    // running list of raw GenomeLocs

        if (argList != null) { // now that we can be in this function if only the ROD-to-Intervals was provided, we need to
            // ensure that the arg list isn't null before looping.
            for (String argument : argList) {
                rawIntervals.addAll(parseIntervalArguments(parser, argument));
            }
        }

        return rawIntervals;
    }

    public static List<GenomeLoc> parseIntervalArguments(GenomeLocParser parser, String arg) {
        List<GenomeLoc> rawIntervals = new ArrayList<>();    // running list of raw GenomeLocs

        if ( arg.indexOf(';') != -1 ) {
            throw new UserException.BadArgumentValue("-L " + arg, "The legacy -L \"interval1;interval2\" syntax " +
                    "is no longer supported. Please use one -L argument for each " +
                    "interval or an interval file instead.");
        }

        if (isUnmapped(arg))
            throw new UserException.BadArgumentValue("-L/-XL", arg, "Currently the only way to view unmapped intervals " +
                    "is to perform a traversal of the entire file without specifying any intervals");
            // if it's a file, add items to raw interval list
        else if (isIntervalFile(arg)) {
            try {
                rawIntervals.addAll(intervalFileToList(parser, arg));
            }
            catch ( UserException.MalformedGenomeLoc e ) {
                throw e;
            }
            catch ( Exception e ) {
                throw new UserException.MalformedFile(new File(arg), "Interval file could not be parsed in any supported format.", e);
            }
        }
        // otherwise treat as an interval -> parse and add to raw interval list
        else {
            rawIntervals.add(parser.parseGenomeLoc(arg));
        }

        return rawIntervals;
    }

    /**
     * Read a file of genome locations to process. The file may be in BED, Picard,
     * or GATK interval format.
     *
     * @param glParser   GenomeLocParser
     * @param file_name  interval file
     * @return List<GenomeLoc> List of Genome Locs that have been parsed from file
     */
    public static List<GenomeLoc> intervalFileToList(final GenomeLocParser glParser, final String file_name) {
        // try to open file
        File inputFile = new File(file_name);
        List<GenomeLoc> ret = new ArrayList<>();

        // case: BED file
        if ( file_name.toUpperCase().endsWith(".BED") ) {
            // this is now supported in Tribble
            throw new UserException("BED files must be parsed through Tribble; parsing them as intervals through the GATK engine is no longer supported");
        }
        else {
            /**
             * IF not a BED file:
             * first try to read it as a Picard interval file since that's well structured
             * we'll fail quickly if it's not a valid file.
             */
            boolean isPicardInterval = false;
            try {
                // Note: Picard will skip over intervals with contigs not in the sequence dictionary
                IntervalList il = IntervalList.fromFile(inputFile);
                isPicardInterval = true;

                int nInvalidIntervals = 0;
                for (Interval interval : il.getIntervals()) {
                    if ( glParser.isValidGenomeLoc(interval.getContig(), interval.getStart(), interval.getEnd(), true))
                        ret.add(glParser.createGenomeLoc(interval.getContig(), interval.getStart(), interval.getEnd(), true));
                    else {
                        nInvalidIntervals++;
                    }
                }
                if ( nInvalidIntervals > 0 )
                    logger.warn("Ignoring " + nInvalidIntervals + " invalid intervals from " + inputFile);
            }

            // if that didn't work, try parsing file as a GATK interval file
            catch (Exception e) {
                if ( isPicardInterval ) // definitely a picard file, but we failed to parse
                    throw new UserException.CouldNotReadInputFile(inputFile, e);
                else {
                    try {
                        XReadLines reader = new XReadLines(new File(file_name));
                        for(String line: reader) {
                            if ( line.trim().length() > 0 ) {
                                ret.add(glParser.parseGenomeLoc(line));
                            }
                        }
                        reader.close();
                    }
                    catch (IOException e2) {
                        throw new UserException.CouldNotReadInputFile(inputFile, e2);
                    }
                }
            }
        }

        if ( ret.isEmpty() ) {
            throw new UserException.MalformedFile(new File(file_name), "It contains no intervals.");
        }

        return ret;
    }

    /**
     * Returns true if the interval string is the "unmapped" interval
     * @param interval Interval to check
     * @return true if the interval string is the "unmapped" interval
     */
    public static boolean isUnmapped(String interval) {
        return (interval != null && interval.trim().toLowerCase().equals("unmapped"));
    }

    /**
     * merge two interval lists, using an interval set rule
     * @param setOne a list of genomeLocs, in order (cannot be NULL)
     * @param setTwo a list of genomeLocs, also in order (cannot be NULL)
     * @param rule the rule to use for merging, i.e. union, intersection, etc
     * @return a list, correctly merged using the specified rule
     */
    public static List<GenomeLoc> mergeListsBySetOperator(List<GenomeLoc> setOne, List<GenomeLoc> setTwo, IntervalSetRule rule) {
        // shortcut, if either set is zero, return the other set
        if (setOne == null || setOne.size() == 0 || setTwo == null || setTwo.size() == 0)
            return Collections.unmodifiableList((setOne == null || setOne.size() == 0) ? setTwo : setOne);

        // our master list, since we can't guarantee removal time in a generic list
        LinkedList<GenomeLoc> retList = new LinkedList<>();

        // if we're set to UNION, just add them all
        if (rule == null || rule == IntervalSetRule.UNION) {
            retList.addAll(setOne);
            retList.addAll(setTwo);
            return Collections.unmodifiableList(retList);
        }

        // else we're INTERSECTION, create two indexes into the lists
        int iOne = 0;
        int iTwo = 0;

        // merge the second into the first using the rule
        while (iTwo < setTwo.size() && iOne < setOne.size())
            // if the first list is ahead, drop items off the second until we overlap
            if (setTwo.get(iTwo).isBefore(setOne.get(iOne)))
                iTwo++;
                // if the second is ahead, drop intervals off the first until we overlap
            else if (setOne.get(iOne).isBefore(setTwo.get(iTwo)))
                iOne++;
                // we overlap, intersect the two intervals and add the result.  Then remove the interval that ends first.
            else {
                retList.add(setOne.get(iOne).intersect(setTwo.get(iTwo)));
                if (setOne.get(iOne).getStop() < setTwo.get(iTwo).getStop()) iOne++;
                else iTwo++;
            }

        //if we have an empty list, throw an exception.  If they specified intersection and there are no items, this is bad.
        if (retList.size() == 0)
            throw new UserException.EmptyIntersection("There was an empty intersection");

        // we don't need to add the rest of remaining locations, since we know they don't overlap. return what we have
        return Collections.unmodifiableList(retList);
    }

    /**
     * Sorts and merges an interval list.  Multiple techniques are available for merging: ALL, which combines
     * all overlapping and abutting intervals into an interval that spans the union of all covered bases, and
     * OVERLAPPING_ONLY, which unions overlapping intervals but keeps abutting intervals separate.
     *
     * @param parser Genome loc parser for the intervals.
     * @param intervals A collection of intervals to merge.
     * @param mergingRule A descriptor for the type of merging to perform.
     * @return A sorted, merged version of the intervals passed in.
     */
    public static GenomeLocSortedSet sortAndMergeIntervals(GenomeLocParser parser, List<GenomeLoc> intervals, IntervalMergingRule mergingRule) {
        // Make a copy of the (potentially unmodifiable) list to be sorted
        intervals = new ArrayList<>(intervals);
        // sort raw interval list
        Collections.sort(intervals);
        // now merge raw interval list
        intervals = mergeIntervalLocations(intervals, mergingRule);

        return GenomeLocSortedSet.createSetFromList(parser, intervals);
    }

    /**
     * computes whether the test interval list is equivalent to master.  To be equivalent, test must
     * contain GenomeLocs covering every base in master, exactly once.  Note that this algorithm
     * assumes that master genomelocs are all discontiguous (i.e., we don't have locs like 1-3 and 4-6 but
     * rather just 1-6).  In order to use this algorithm with contiguous genomelocs first merge them.  The algorithm
     * doesn't assume that test has discontinuous genomelocs.
     *
     * Returns a null string if there are no differences, otherwise returns a string describing the difference
     * (useful for UnitTests).  Assumes both lists are sorted
     *
     * @param masterArg sorted master genome locs
     * @param testArg sorted test genome locs
     * @return null string if there are no difference, otherwise a string describing the difference
     */
    public static String equateIntervals(List<GenomeLoc> masterArg, List<GenomeLoc> testArg) {
        LinkedList<GenomeLoc> master = new LinkedList<>(masterArg);
        LinkedList<GenomeLoc> test = new LinkedList<>(testArg);

        while ( ! master.isEmpty() ) { // there's still unchecked bases in master
            final GenomeLoc masterHead = master.pop();
            final GenomeLoc testHead = test.pop();

            if ( testHead.overlapsP(masterHead) ) {
                // remove the parts of test that overlap master, and push the remaining
                // parts onto master for further comparison.
                Utils.reverse(masterHead.subtract(testHead)).forEach(master::push);
            } else {
                // testHead is incompatible with masterHead, so we must have extra bases in testHead
                // that aren't in master
                return "Incompatible locs detected masterHead=" + masterHead + ", testHead=" + testHead;
            }
        }

        if ( test.isEmpty() ) // everything is equal
            return null; // no differences
        else
            return "Remaining elements found in test: first=" + test.peek();
    }


    /**
     * Check if string argument was intented as a file
     * Accepted file extensions: .bed .list, .picard, .interval_list, .intervals.
     * @param str token to identify as a filename.
     * @return true if the token looks like a filename, or false otherwise.
     */
    public static boolean isIntervalFile(String str) {
        return isIntervalFile(str, true);
    }

    /**
     * Check if string argument was intented as a file
     * Accepted file extensions: .bed .list, .picard, .interval_list, .intervals.
     * @param str token to identify as a filename.
     * @param checkExists if true throws an exception if the file doesn't exist.
     * @return true if the token looks like a filename, or false otherwise.
     */
    public static boolean isIntervalFile(String str, boolean checkExists) {
        // should we define list of file extensions as a public array somewhere?
        // is regex or endsiwth better?
        File file = new File(str);
        if (str.toUpperCase().endsWith(".BED") || str.toUpperCase().endsWith(".LIST") ||
                str.toUpperCase().endsWith(".PICARD") || str.toUpperCase().endsWith(".INTERVAL_LIST")
                || str.toUpperCase().endsWith(".INTERVALS")) {
            if (!checkExists)
                return true;
            else if (file.exists())
                return true;
            else
                throw new UserException.CouldNotReadInputFile(file, "The interval file does not exist.");
        }

        if(file.exists())
            throw new UserException.CouldNotReadInputFile(file, String.format("The interval file %s does not have one of " +
                    "the supported extensions (.bed, .list, .picard, .interval_list, or .intervals). " +
                    "Please rename your file with the appropriate extension. If %s is NOT supposed to be a file, " +
                    "please move or rename the file at location %s", str, str, file.getAbsolutePath()));

        else return false;
    }

    /**
     * Returns a map of contig names with their sizes.
     * @param reference The reference for the intervals.
     * @return A map of contig names with their sizes.
     */
    public static Map<String, Integer> getContigSizes(File reference) {
        final ReferenceSequenceFile referenceSequenceFile = createReference(reference);
        List<GenomeLoc> locs = GenomeLocSortedSet.createSetFromSequenceDictionary(referenceSequenceFile.getSequenceDictionary()).toList();
        Map<String, Integer> lengths = new LinkedHashMap<>();
        for (GenomeLoc loc: locs)
            lengths.put(loc.getContig(), loc.size());
        return lengths;
    }

    /**
     * Splits an interval list into multiple files.
     * @param fileHeader The sam file header.
     * @param locs The genome locs to split.
     * @param scatterParts The output interval lists to write to.
     */
    public static void scatterContigIntervals(SAMFileHeader fileHeader, List<GenomeLoc> locs, List<File> scatterParts) {

        // Contract: must divide locs up so that each of scatterParts gets a sublist such that:
        // (a) all locs concerning a particular contig go to the same part
        // (b) locs are not split or combined, and remain in the same order (so scatterParts[0] + ... + scatterParts[n] == locs)

        // Locs are already sorted.

        long totalBases = 0;
        for(GenomeLoc loc : locs)
            totalBases += loc.size();

        long idealBasesPerPart = totalBases / scatterParts.size();
        if(idealBasesPerPart == 0)
            throw new UserException.BadInput(String.format("Genome region is too short (%d bases) to split into %d parts", totalBases, scatterParts.size()));

        // Find the indices in locs where we switch from one contig to the next.
        ArrayList<Integer> contigStartLocs = new ArrayList<>();
        String prevContig = null;

        for(int i = 0; i < locs.size(); ++i) {

            GenomeLoc loc = locs.get(i);
            if(prevContig == null || !loc.getContig().equals(prevContig))
                contigStartLocs.add(i);
            prevContig = loc.getContig();

        }

        if(contigStartLocs.size() < scatterParts.size())
            throw new UserException.BadInput(String.format("Input genome region has too few contigs (%d) to split into %d parts", contigStartLocs.size(), scatterParts.size()));

        long thisPartBases = 0;
        int partIdx = 0;
        IntervalList outList = new IntervalList(fileHeader);

        for(int i = 0; i < locs.size(); ++i) {

            GenomeLoc loc = locs.get(i);
            thisPartBases += loc.getStop() - loc.getStart();

            outList.add(toInterval(loc, i));

            boolean partMustStop = false;

            if(partIdx < (scatterParts.size() - 1)) {

                // If there are n contigs and n parts remaining then we must split here,
                // otherwise we will run out of contigs.

                int nextPart = partIdx + 1;
                int nextPartMustStartBy = contigStartLocs.get(nextPart + (contigStartLocs.size() - scatterParts.size()));
                if(i + 1 == nextPartMustStartBy)
                    partMustStop = true;

            }
            else if(i == locs.size() - 1) {

                // We're done! Write the last scatter file.
                partMustStop = true;

            }

            if(partMustStop || thisPartBases > idealBasesPerPart) {

                // Ideally we would split here. However, we must make sure to do so
                // on a contig boundary. Test always passes with partMustStop == true
                // since that indicates we're at a contig boundary.

                GenomeLoc nextLoc = null;
                if((i + 1) < locs.size())
                    nextLoc = locs.get(i+1);

                if(nextLoc == null || !nextLoc.getContig().equals(loc.getContig())) {

                    // Write out this part:
                    outList.write(scatterParts.get(partIdx));

                    // Reset. If this part ran long, leave the excess in thisPartBases
                    // and the next will be a little shorter to compensate.
                    outList = new IntervalList(fileHeader);
                    thisPartBases -= idealBasesPerPart;
                    ++partIdx;

                }

            }

        }

    }

    /**
     * Splits an interval list into multiple sublists.
     * @param locs The genome locs to split.
     * @param splits The stop points for the genome locs returned by splitFixedIntervals.
     * @return A list of lists of genome locs, split according to splits
     */
    public static List<List<GenomeLoc>> splitIntervalsToSubLists(List<GenomeLoc> locs, List<Integer> splits) {
        int start = 0;
        List<List<GenomeLoc>> sublists = new ArrayList<>(splits.size());
        for (Integer stop: splits) {
            List<GenomeLoc> curList = new ArrayList<>();
            for (int i = start; i < stop; i++)
                curList.add(locs.get(i));
            start = stop;
            sublists.add(curList);
        }

        return sublists;
    }


    /**
     * Splits an interval list into multiple files.
     * @param fileHeader The sam file header.
     * @param splits Pre-divided genome locs returned by splitFixedIntervals.
     * @param scatterParts The output interval lists to write to.
     */
    public static void scatterFixedIntervals(SAMFileHeader fileHeader, List<List<GenomeLoc>> splits, List<File> scatterParts) {
        if (splits.size() != scatterParts.size())
            throw new UserException.BadArgumentValue("splits", String.format("Split points %d does not equal the number of scatter parts %d.", splits.size(), scatterParts.size()));

        int fileIndex = 0;
        int locIndex = 1;
        for (final List<GenomeLoc> split : splits) {
            IntervalList intervalList = new IntervalList(fileHeader);
            for (final GenomeLoc loc : split)
                intervalList.add(toInterval(loc, locIndex++));
            intervalList.write(scatterParts.get(fileIndex++));
        }
    }

    /**
     * Splits the genome locs up by size.
     * @param locs Genome locs to split.
     * @param numParts Number of parts to split the locs into.
     * @return The stop points to split the genome locs.
     */
    public static List<List<GenomeLoc>> splitFixedIntervals(List<GenomeLoc> locs, int numParts) {
        if (locs.size() < numParts)
            throw new UserException.BadArgumentValue("scatterParts", String.format("Cannot scatter %d locs into %d parts.", locs.size(), numParts));
        final long locsSize = intervalSize(locs);
        final List<Integer> splitPoints = new ArrayList<>();
        addFixedSplit(splitPoints, locs, locsSize, 0, locs.size(), numParts);
        Collections.sort(splitPoints);
        splitPoints.add(locs.size());
        return splitIntervalsToSubLists(locs, splitPoints);
    }

    public static List<List<GenomeLoc>> splitLocusIntervals(List<GenomeLoc> locs, int numParts) {
        // the ideal size of each split
        final long bp = IntervalUtils.intervalSize(locs);
        final long idealSplitSize = Math.max((long)Math.floor(bp / (1.0*numParts)), 1);

        // algorithm:
        // split = ()
        // set size = 0
        // pop the head H off locs.
        // If size + size(H) < splitSize:
        //      add H to split, continue
        // If size + size(H) == splitSize:
        //      done with split, put in splits, restart
        // if size + size(H) > splitSize:
        //      cut H into two pieces, first of which has splitSize - size bp
        //      push both pieces onto locs, continue
        // The last split is special -- when you have only one split left, it gets all of the remaining locs
        // to deal with rounding issues
        final List<List<GenomeLoc>> splits = new ArrayList<>(numParts);

        LinkedList<GenomeLoc> locsLinkedList = new LinkedList<>(locs);
        while ( ! locsLinkedList.isEmpty() ) {
            if ( splits.size() + 1 == numParts ) {
                // the last one gets all of the remaining parts
                splits.add(new ArrayList<>(locsLinkedList));
                locsLinkedList.clear();
            } else {
                final SplitLocusRecursive one = splitLocusIntervals1(locsLinkedList, idealSplitSize);
                splits.add(one.split);
                locsLinkedList = one.remaining;
            }
        }

        return splits;
    }

    static SplitLocusRecursive splitLocusIntervals1(LinkedList<GenomeLoc> remaining, long idealSplitSize) {
        final List<GenomeLoc> split = new ArrayList<>();
        long size = 0;

        while ( ! remaining.isEmpty() ) {
            GenomeLoc head = remaining.pop();
            final long newSize = size + head.size();

            if ( newSize == idealSplitSize ) {
                split.add(head);
                break; // we are done
            } else if ( newSize > idealSplitSize ) {
                final long remainingBp = idealSplitSize - size;
                final long cutPoint = head.getStart() + remainingBp;
                GenomeLoc[] parts = head.split((int)cutPoint);
                remaining.push(parts[1]);
                remaining.push(parts[0]);
                // when we go around, head.size' = idealSplitSize - size
                // so newSize' = splitSize + head.size' = size + (idealSplitSize - size) = idealSplitSize
            } else {
                split.add(head);
                size = newSize;
            }
        }

        return new SplitLocusRecursive(split, remaining);
    }

    /**
     * Check whether two locatables overlap.
     * <p>
     *    Two locatables overlap if the share the same contig and they have at least one
     *    base in common based on their start and end positions.
     * </p>
     * <p>
     *    This method returns {@code false} if either input {@link Locatable} has a {@code null}
     *    contig.
     * </p>
     *
     * @param left first locatable.
     * @param right second locatable.
     * @throws IllegalArgumentException if either {@code left} or {@code right} locatable
     *  is {@code null}.
     * @return {@code true} iff there is an overlap between both locatables.
     */
    public static boolean overlaps(final Locatable left,final Locatable right) {
        Utils.nonNull(left,"the left locatable is null");
        Utils.nonNull(right,"the right locatable is null");
        if (left.getContig() == null || right.getContig() == null) {
            return false;
        } else if (!left.getContig().equals(right.getContig())) {
            return false;
        } else {
            return left.getStart() <= right.getEnd() && right.getStart() <= left.getEnd();
        }
    }

    private final static class SplitLocusRecursive {
        final List<GenomeLoc> split;
        final LinkedList<GenomeLoc> remaining;

        private SplitLocusRecursive(final List<GenomeLoc> split, final LinkedList<GenomeLoc> remaining) {
            this.split = split;
            this.remaining = remaining;
        }
    }

    public static List<GenomeLoc> flattenSplitIntervals(List<List<GenomeLoc>> splits) {
        final List<GenomeLoc> locs = new ArrayList<>();
        splits.forEach(locs::addAll);

        return locs;
    }

    private static void addFixedSplit(List<Integer> splitPoints, List<GenomeLoc> locs, long locsSize, int startIndex, int stopIndex, int numParts) {
        if (numParts < 2)
            return;
        int halfParts = (numParts + 1) / 2;
        Pair<Integer, Long> splitPoint = getFixedSplit(locs, locsSize, startIndex, stopIndex, halfParts, numParts - halfParts);
        int splitIndex = splitPoint.getLeft();
        long splitSize = splitPoint.getRight();
        splitPoints.add(splitIndex);
        addFixedSplit(splitPoints, locs, splitSize, startIndex, splitIndex, halfParts);
        addFixedSplit(splitPoints, locs, locsSize - splitSize, splitIndex, stopIndex, numParts - halfParts);
    }

    private static Pair<Integer, Long> getFixedSplit(List<GenomeLoc> locs, long locsSize, int startIndex, int stopIndex, int minLocs, int maxLocs) {
        int splitIndex = startIndex;
        long splitSize = 0;
        for (int i = 0; i < minLocs; i++) {
            splitSize += locs.get(splitIndex).size();
            splitIndex++;
        }
        long halfSize = locsSize / 2;
        while (splitIndex < (stopIndex - maxLocs) && splitSize < halfSize) {
            splitSize += locs.get(splitIndex).size();
            splitIndex++;
        }
        return new ImmutablePair<>(splitIndex, splitSize);
    }

    /**
     * Converts a GenomeLoc to a picard interval.
     * @param loc The GenomeLoc.
     * @param locIndex The loc index for use in the file.
     * @return The picard interval.
     */
    private static htsjdk.samtools.util.Interval toInterval(GenomeLoc loc, int locIndex) {
        return new htsjdk.samtools.util.Interval(loc.getContig(), loc.getStart(), loc.getStop(), false, "interval_" + locIndex);
    }

    /**
     * merge a list of genome locs that may be overlapping, returning the list of unique genomic locations
     *
     * @param raw the unchecked genome loc list
     * @param rule the merging rule we're using
     *
     * @return the list of merged locations
     */
    public static List<GenomeLoc> mergeIntervalLocations(final List<GenomeLoc> raw, IntervalMergingRule rule) {
        if (raw.size() <= 1)
            return Collections.unmodifiableList(raw);
        else {
            ArrayList<GenomeLoc> merged = new ArrayList<>();
            Iterator<GenomeLoc> it = raw.iterator();
            GenomeLoc prev = it.next();
            while (it.hasNext()) {
                GenomeLoc curr = it.next();
                if (prev.overlapsP(curr)) {
                    prev = prev.merge(curr);
                } else if (prev.contiguousP(curr) && (rule == null || rule == IntervalMergingRule.ALL)) {
                    prev = prev.merge(curr);
                } else {
                    merged.add(prev);
                    prev = curr;
                }
            }
            merged.add(prev);
            return Collections.unmodifiableList(merged);
        }
    }

    public static long intervalSize(final List<GenomeLoc> locs) {
        long size = 0;
        for ( final GenomeLoc loc : locs )
            size += loc.size();
        return size;
    }

    public static void writeFlankingIntervals(File reference, File inputIntervals, File flankingIntervals, int basePairs) {
        final ReferenceSequenceFile referenceSequenceFile = createReference(reference);
        GenomeLocParser parser = new GenomeLocParser(referenceSequenceFile);
        List<GenomeLoc> originalList = intervalFileToList(parser, inputIntervals.getAbsolutePath());

        if (originalList.isEmpty())
            throw new UserException.MalformedFile(inputIntervals, "File contains no intervals");

        List<GenomeLoc> flankingList = getFlankingIntervals(parser, originalList, basePairs);

        if (flankingList.isEmpty())
            throw new UserException.MalformedFile(inputIntervals, "Unable to produce any flanks for the intervals");

        SAMFileHeader samFileHeader = new SAMFileHeader();
        samFileHeader.setSequenceDictionary(referenceSequenceFile.getSequenceDictionary());
        IntervalList intervalList = new IntervalList(samFileHeader);
        int i = 0;
        for (GenomeLoc loc: flankingList)
            intervalList.add(toInterval(loc, ++i));
        intervalList.write(flankingIntervals);
    }

    /**
     * Returns a list of intervals between the passed int locs. Does not extend UNMAPPED locs.
     * @param parser A genome loc parser for creating the new intervals
     * @param locs Original genome locs
     * @param basePairs Number of base pairs on each side of loc
     * @return The list of intervals between the locs
     */
    public static List<GenomeLoc> getFlankingIntervals(final GenomeLocParser parser, final List<GenomeLoc> locs, final int basePairs) {
        List<GenomeLoc> sorted = sortAndMergeIntervals(parser, locs, IntervalMergingRule.ALL).toList();

        if (sorted.size() == 0)
            return Collections.emptyList();

        LinkedHashMap<String, List<GenomeLoc>> locsByContig = splitByContig(sorted);
        List<GenomeLoc> expanded = new ArrayList<>();
        for (Map.Entry<String, List<GenomeLoc>> contig: locsByContig.entrySet()) {
            List<GenomeLoc> contigLocs = contig.getValue();
            int contigLocsSize = contigLocs.size();

            GenomeLoc startLoc, stopLoc;

            // Create loc at start of the list
            startLoc = parser.createGenomeLocAtStart(contigLocs.get(0), basePairs);
            if (startLoc != null)
                expanded.add(startLoc);

            // Create locs between each loc[i] and loc[i+1]
            for (int i = 0; i < contigLocsSize - 1; i++) {
                stopLoc = parser.createGenomeLocAtStop(contigLocs.get(i), basePairs);
                startLoc = parser.createGenomeLocAtStart(contigLocs.get(i + 1), basePairs);
                if (stopLoc.getStop() + 1 >= startLoc.getStart()) {
                    // NOTE: This is different than GenomeLoc.merge()
                    // merge() returns a loc which covers the entire range of stop and start,
                    // possibly returning positions inside loc(i) or loc(i+1)
                    // We want to make sure that the start of the stopLoc is used, and the stop of the startLoc
                    GenomeLoc merged = parser.createGenomeLoc(
                            stopLoc.getContig(), stopLoc.getStart(), startLoc.getStop());
                    expanded.add(merged);
                } else {
                    expanded.add(stopLoc);
                    expanded.add(startLoc);
                }
            }

            // Create loc at the end of the list
            stopLoc = parser.createGenomeLocAtStop(contigLocs.get(contigLocsSize - 1), basePairs);
            if (stopLoc != null)
                expanded.add(stopLoc);
        }
        return expanded;
    }

    /**
     * Returns a list of intervals between the passed int locs. Does not extend UNMAPPED locs.
     * @param parser A genome loc parser for creating the new intervals
     * @param locs Original genome locs
     * @param basePairs Number of base pairs on each side of loc
     * @return The list of intervals between the locs
     */
    public static List<GenomeLoc> getIntervalsWithFlanks(final GenomeLocParser parser, final List<GenomeLoc> locs, final int basePairs) {

        if (locs.size() == 0)
            return Collections.emptyList();

        final List<GenomeLoc> expanded = locs.stream()
                .map(loc -> parser.createPaddedGenomeLoc(loc, basePairs))
                .collect(Collectors.toList());

        return sortAndMergeIntervals(parser, expanded, IntervalMergingRule.ALL).toList();
    }

    private static ReferenceSequenceFile createReference(final File fastaFile) {
            return CachingIndexedFastaSequenceFile.checkAndCreate(fastaFile);
    }

    private static LinkedHashMap<String, List<GenomeLoc>> splitByContig(List<GenomeLoc> sorted) {
        LinkedHashMap<String, List<GenomeLoc>> splits = new LinkedHashMap<>();
        GenomeLoc last = null;
        List<GenomeLoc> contigLocs = null;
        for (GenomeLoc loc: sorted) {
            if (GenomeLoc.isUnmapped(loc))
                continue;
            if (last == null || !last.onSameContig(loc)) {
                contigLocs = new ArrayList<>();
                splits.put(loc.getContig(), contigLocs);
            }
            contigLocs.add(loc);
            last = loc;
        }
        return splits;
    }

    /**
     * Generates a list of {@link GenomeLoc} instances given the appropriate {@link GenomeLocParser} factory
     * and a collection of {@link Locatable} instances.
     *
     * <p>
     *     The order in the result list is will correspond to the traversal order in the input collection.
     * </p>
     * @param locatables input locatable collection.
     * @throws IllegalArgumentException if {@code locatable} is {@code null} or contains any {@code null}.
     * @return never {@code null}. The result is an unmodifiable list.
     */
    public static List<GenomeLoc> genomeLocsFromLocatables(final GenomeLocParser parser, final Collection<? extends Locatable> locatables) {
        if (parser == null) {
            throw new IllegalArgumentException("the input genome-loc parser cannot be null");
        }
        if (locatables == null) {
            throw new IllegalArgumentException("the input locatable collection cannot be null");
        };
        if (locatables.stream().anyMatch(Objects::isNull)) {
            throw new IllegalArgumentException("some element in the locatable input collection is null");
        }
        final List<GenomeLoc> result = locatables.stream().map(parser::createGenomeLoc).collect(Collectors.toList());
        return Collections.unmodifiableList(result);
    }

    /**
     * Builds a list of intervals that cover the whole given sequence.
     */
    public static List<SimpleInterval> getAllIntervalsForReference(SAMSequenceDictionary sequenceDictionary) {
        return GenomeLocSortedSet.createSetFromSequenceDictionary(sequenceDictionary)
                .stream()
                .map(SimpleInterval::new)
                .collect(Collectors.toList());
    }

}
