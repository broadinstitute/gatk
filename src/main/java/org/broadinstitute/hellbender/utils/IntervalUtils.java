package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.Feature;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.nio.PathLineIterator;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Parse text representations of interval strings that
 * can appear in GATK-based applications.
 */
public final class IntervalUtils {

    private static final Logger log = LogManager.getLogger(IntervalUtils.class);

    /**
     * Recognized extensions for interval files
     */
    public static final List<String> INTERVAL_FILE_EXTENSIONS = Collections.unmodifiableList(Arrays.asList(
        ".list", ".interval_list", ".intervals", ".picard"
    ));

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
    public static final Comparator<Locatable> LEXICOGRAPHICAL_ORDER_COMPARATOR =
            Comparator.comparing(Locatable::getContig,Comparator.nullsLast(String::compareTo))
                    .thenComparingInt(Locatable::getStart)
                    .thenComparingInt(Locatable::getEnd);

    private static final Logger logger = LogManager.getLogger(IntervalUtils.class);

    /**
     * Compare two locations using a {@link htsjdk.samtools.SAMSequenceDictionary} sequence ordering
     *
     * @return 0 if first and second are equal, a negative value if first < second or a positive value if first > second,
     * using the ordering in the provided dictionary
     * @throws IllegalArgumentException if either first or second contigs could not be found in the dictionary
     */
    public static final int compareLocatables(Locatable first, Locatable second, SAMSequenceDictionary dictionary) {
        Utils.nonNull(first);
        Utils.nonNull(second);
        Utils.nonNull(dictionary);

        int result = 0;
        if(first != second) {
            // compare the contigs
            result = compareContigs(first, second, dictionary);
            if (result == 0) {
                // compare start position
                result = Integer.compare(first.getStart(), second.getStart());
                if (result == 0) {
                    // compare end position
                    result = Integer.compare(first.getEnd(), second.getEnd());
                }
            }
        }
        return result;
    }

    /**
     * Tests whether the first Locatable ends before the start of the second Locatable
     *
     * @param first first Locatable
     * @param second second Locatable
     * @param dictionary sequence dictionary used to determine contig ordering
     * @return true if first ends before the start of second, otherwise false
     */
    public static boolean isBefore(final Locatable first, final Locatable second, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(first);
        Utils.nonNull(second);
        Utils.nonNull(dictionary);

        final int contigComparison = compareContigs(first, second, dictionary);
        return contigComparison == -1 || (contigComparison == 0 && first.getEnd() < second.getStart());
    }

    /**
     * Tests whether the first Locatable starts after the end of the second Locatable
     *
     * @param first first Locatable
     * @param second second Locatable
     * @param dictionary sequence dictionary used to determine contig ordering
     * @return true if first starts after the end of second, otherwise false
     */
    public static boolean isAfter(final Locatable first, final Locatable second, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(first);
        Utils.nonNull(second);
        Utils.nonNull(dictionary);

        final int contigComparison = compareContigs(first, second, dictionary);
        return contigComparison == 1 || (contigComparison == 0 && first.getStart() > second.getEnd());
    }

    /**
     * Determines the relative contig ordering of first and second using the provided sequence dictionary
     *
     * @param first first Locatable
     * @param second second Locatable
     * @param dictionary sequence dictionary used to determine contig ordering
     * @return 0 if the two contigs are the same, a negative value if first's contig comes before second's contig,
     *         or a positive value if first's contig comes after second's contig
     */
    public static int compareContigs(final Locatable first, final Locatable second, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(first);
        Utils.nonNull(second);
        Utils.nonNull(dictionary);

        final int firstRefIndex = dictionary.getSequenceIndex(first.getContig());
        final int secondRefIndex = dictionary.getSequenceIndex(second.getContig());
        if (firstRefIndex == -1 || secondRefIndex == -1) {
            throw new IllegalArgumentException("Can't do comparison because Locatables' contigs not found in sequence dictionary");
        }
        // compare the contigs
        return Integer.compare(firstRefIndex, secondRefIndex);
    }


    /**
     * getSpanningInterval returns interval that covers all of the locations passed in.
     * @param locations the locations to be spanned (on a single contig)
     * @return the minimal span that covers all locations (could be null if no locations are passed in).
     * @throws IllegalArgumentException if the argument is null
     * or if the argument contains any null element
     * or if the locations are not all on the same contig (compared by String.equals)
     */
    public static SimpleInterval getSpanningInterval(final List<? extends Locatable> locations) {
        Utils.nonNull(locations);
        Utils.containsNoNull(locations, "locations must not contain a null");
        if (locations.isEmpty()){
            return null;
        }
        final List<String> contigs = locations.stream().map(l -> l.getContig()).distinct().collect(Collectors.toList());
        Utils.validateArg(contigs.size() == 1, () -> "found different contigs from inputs:" + contigs);
        final int minStart = locations.stream().mapToInt(l -> l.getStart()).min().getAsInt();
        final int maxEnd   = locations.stream().mapToInt(l -> l.getEnd()).max().getAsInt();
        return new SimpleInterval(contigs.get(0), minStart, maxEnd);
    }

    /**
     * Convert a List of intervals in GenomeLoc format into a List of intervals in SimpleInterval format.
     *
     * @param genomeLocIntervals list of GenomeLoc intervals to convert
     * @return equivalent List of SimpleIntervals
     */
    public static List<SimpleInterval> convertGenomeLocsToSimpleIntervals( final List<GenomeLoc> genomeLocIntervals ) {
        final List<SimpleInterval> convertedIntervals = new ArrayList<>(genomeLocIntervals.size());
        for ( final GenomeLoc genomeLoc : genomeLocIntervals ) {
            if ( genomeLoc.isUnmapped() ) {
                throw new UserException("Unmapped intervals cannot be converted to SimpleIntervals");
            }

            convertedIntervals.add(new SimpleInterval(genomeLoc));
        }
        return convertedIntervals;
    }

    /**
     * Converts an interval in SimpleInterval format into an htsjdk QueryInterval.
     *
     * In doing so, a header lookup is performed to convert from contig name to index
     *
     * @param interval interval to convert
     * @param sequenceDictionary sequence dictionary used to perform the conversion
     * @return an equivalent interval in QueryInterval format
     */
    public static QueryInterval convertSimpleIntervalToQueryInterval( final SimpleInterval interval, final SAMSequenceDictionary sequenceDictionary ) {
        Utils.nonNull(interval);
        Utils.nonNull(sequenceDictionary);

        final int contigIndex = sequenceDictionary.getSequenceIndex(interval.getContig());
        if ( contigIndex == -1 ) {
            throw new UserException("Contig " + interval.getContig() + " not present in reads sequence dictionary");
        }

        return new QueryInterval(contigIndex, interval.getStart(), interval.getEnd());
    }

    public static GenomeLocSortedSet loadIntervals(
            final List<String> intervalStrings,
            final IntervalSetRule intervalSetRule,
            final IntervalMergingRule intervalMergingRule,
            final int padding,
            final GenomeLocParser genomeLocParser) {
        Utils.nonNull(intervalStrings);
        List<GenomeLoc> allIntervals = new ArrayList<>();
        for ( final String intervalString : intervalStrings) {
            Utils.nonNull(intervalString);
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
    public static List<GenomeLoc> parseIntervalArguments(final GenomeLocParser parser, final List<String> argList) {
        final List<GenomeLoc> rawIntervals = new ArrayList<>();    // running list of raw GenomeLocs

        if (argList != null) { // now that we can be in this function if only the ROD-to-Intervals was provided, we need to
            // ensure that the arg list isn't null before looping.
            for (final String argument : argList) {
                rawIntervals.addAll(parseIntervalArguments(parser, argument));
            }
        }

        return rawIntervals;
    }

    public static List<GenomeLoc> parseIntervalArguments(final GenomeLocParser parser, final String arg) {
        Utils.nonNull(parser, "parser is null");
        Utils.nonNull(arg, "arg is null");
        final List<GenomeLoc> rawIntervals = new ArrayList<>();    // running list of raw GenomeLocs

        if ( arg.indexOf(';') != -1 ) {
            throw new CommandLineException.BadArgumentValue("-L " + arg, "The legacy -L \"interval1;interval2\" syntax " +
                    "is no longer supported. Please use one -L argument for each " +
                    "interval or an interval file instead.");
        }
        // If it's a Feature-containing file, convert it to a list of intervals
        else if ( FeatureManager.isFeatureFile(IOUtils.getPath(arg)) ) {
            rawIntervals.addAll(featureFileToIntervals(parser, arg));
        }
        // If it's an interval file, add its contents to the raw interval list
        else if ( isIntervalFile(arg) ) {
            try {
                rawIntervals.addAll(intervalFileToList(parser, arg));
            }
            catch ( final UserException.MalformedGenomeLoc e ) {
                throw e;
            }
            catch ( final Exception e ) {
                throw new UserException.MalformedFile(new File(arg), "Interval file could not be parsed in any supported format.", e);
            }
        }
        // If it's neither a Feature-containing file nor an interval file, but is an existing file, throw an error.
        // Note that since contigs can contain periods in their names, we can't use the mere presence of an "extension"
        // as evidence that the user intended the String to be interpreted as a file.
        else if ( new File(arg).exists() ) {
            throw new UserException.CouldNotReadInputFile(arg, String.format("The file %s exists, but does not contain Features " +
                    "(ie., is not in a supported Feature file format such as vcf, bcf, or bed), " +
                    "and does not have one of the supported interval file extensions (" + INTERVAL_FILE_EXTENSIONS + "). " +
                    "Please rename your file with the appropriate extension. If %s is NOT supposed to be a file, " +
                    "please move or rename the file at location %s", arg, arg, new File(arg).getAbsolutePath()));
        }
        // Otherwise treat as an interval -> parse and add to raw interval list
        else {
            rawIntervals.add(parser.parseGenomeLoc(arg));
        }

        return rawIntervals;
    }

    /**
     * Converts a Feature-containing file to a list of intervals
     *
     * @param parser GenomeLocParser for creating intervals
     * @param featureFile file containing Features to convert to intervals
     * @return a List of intervals corresponding to the locations of the Features in the provided file
     * @throws UserException.CouldNotReadInputFile if the provided file is not in a supported Feature file format
     */
    public static List<GenomeLoc> featureFileToIntervals( final GenomeLocParser parser, final String featureFile ) {

        try ( final FeatureDataSource<? extends Feature> dataSource = new FeatureDataSource<>(featureFile) ) {
            final List<GenomeLoc> featureIntervals = new ArrayList<>();

            for ( final Feature feature : dataSource ) {
                featureIntervals.add(parser.createGenomeLoc(feature));
            }
            return featureIntervals;
        }
    }

    /**
     * Read a file of genome locations to process. The file may be in Picard
     * or GATK interval format.
     *
     * @param glParser   GenomeLocParser
     * @param fileName  interval file
     * @return List<GenomeLoc> List of Genome Locs that have been parsed from file
     */
    public static List<GenomeLoc> intervalFileToList(final GenomeLocParser glParser, final String fileName) {
        Utils.nonNull(glParser, "glParser is null");
        Utils.nonNull(fileName, "file name is null");

        final Path inputPath = IOUtils.getPath(fileName);
        final List<GenomeLoc> ret = new ArrayList<>();

        /**
         * First try to read the file as a Picard interval file since that's well structured --
         * we'll fail quickly if it's not a valid file.
         */
        boolean isPicardInterval = false;
        try {
            // Note: Picard will skip over intervals with contigs not in the sequence dictionary
            final IntervalList il = IntervalList.fromPath(inputPath);
            isPicardInterval = true;

            for (final Interval interval : il.getIntervals()) {
                // The current Agilent exome interval list is off-by-one on all end positions. Until this is fixed we
                // need to tolerate intervals where the end is one before the start. We should remove this once a
                // corrected version of the interval list is released. This is tracked in:
                // https://github.com/broadinstitute/gatk/issues/2089
                if (interval.getStart() - interval.getEnd() == 1 ) {
                    logger.warn("Ignoring possibly incorrectly converted length 1 interval : " + interval);
                }
                else if ( glParser.isValidGenomeLoc(interval.getContig(), interval.getStart(), interval.getEnd(), true)) {
                    ret.add(glParser.createGenomeLoc(interval.getContig(), interval.getStart(), interval.getEnd(), true));
                } else {
                    throw new UserException(inputPath.toUri() +  " has an invalid interval : " + interval) ;
                }
            }
        }
        // if that didn't work, try parsing file as a GATK interval file
        catch (final Exception e) {
            if ( isPicardInterval ) // definitely a picard file, but we failed to parse
            {
                throw new UserException.CouldNotReadInputFile(inputPath, e);
            } else {
                try (PathLineIterator reader = new PathLineIterator(inputPath)) {
                    for (final String line : reader) {
                        if (!line.trim().isEmpty()) {
                            ret.add(glParser.parseGenomeLoc(line));
                        }
                    }
                }
            }
        }

        if ( ret.isEmpty() ) {
            throw new UserException.MalformedFile(inputPath, "It contains no intervals.");
        }

        return ret;
    }

    /**
     * merge two interval lists, using an interval set rule
     * @param setOne a list of genomeLocs, in order (cannot be NULL)
     * @param setTwo a list of genomeLocs, also in order (cannot be NULL)
     * @param rule the rule to use for merging, i.e. union, intersection, etc
     * @return a list, correctly merged using the specified rule
     */
    public static List<GenomeLoc> mergeListsBySetOperator(final List<GenomeLoc> setOne, final List<GenomeLoc> setTwo, final IntervalSetRule rule) {
        // shortcut, if either set is zero, return the other set
        if (setOne == null || setOne.isEmpty() || setTwo == null || setTwo.isEmpty()) {
            return Collections.unmodifiableList((setOne == null || setOne.isEmpty()) ? setTwo : setOne);
        }

        // our master list, since we can't guarantee removal time in a generic list
        final LinkedList<GenomeLoc> retList = new LinkedList<>();

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
        {
            if (setTwo.get(iTwo).isBefore(setOne.get(iOne))) {
                iTwo++;
            }// if the second is ahead, drop intervals off the first until we overlap
            else if (setOne.get(iOne).isBefore(setTwo.get(iTwo))) {
                iOne++;
            }// we overlap, intersect the two intervals and add the result.  Then remove the interval that ends first.
            else {
                retList.add(setOne.get(iOne).intersect(setTwo.get(iTwo)));
                if (setOne.get(iOne).getStop() < setTwo.get(iTwo).getStop()) {
                    iOne++;
                } else {
                    iTwo++;
                }
            }
        }

        //if we have an empty list, throw an exception.  If they specified intersection and there are no items, this is bad.
        if (retList.isEmpty()) {
            throw new UserException.EmptyIntersection("There was an empty intersection");
        }

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
    public static GenomeLocSortedSet sortAndMergeIntervals(final GenomeLocParser parser, List<GenomeLoc> intervals, final IntervalMergingRule mergingRule) {
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
    public static String equateIntervals(final List<GenomeLoc> masterArg, final List<GenomeLoc> testArg) {
        final LinkedList<GenomeLoc> master = new LinkedList<>(masterArg);
        final LinkedList<GenomeLoc> test = new LinkedList<>(testArg);

        while ( ! master.isEmpty() ) { // there's still unchecked bases in master
            final GenomeLoc masterHead = master.pop();
            final GenomeLoc testHead = test.pop();

            if ( testHead.overlapsP(masterHead) ) {
                // remove the parts of test that overlap master, and push the remaining
                // parts onto master for further comparison.
                reverse(masterHead.subtract(testHead)).forEach(master::push);
            } else {
                // testHead is incompatible with masterHead, so we must have extra bases in testHead
                // that aren't in master
                return "Incompatible locs detected masterHead=" + masterHead + ", testHead=" + testHead;
            }
        }

        if ( test.isEmpty() ) // everything is equal
        {
            return null; // no differences
        } else {
            return "Remaining elements found in test: first=" + test.peek();
        }
    }

    private static <T extends Comparable<T>> List<T> reverse(final List<T> l) {
        return l.stream().sorted(Collections.reverseOrder()).collect(Collectors.toList());
    }

    /**
     * Check if string argument was intended as a file
     * Accepted file extensions: .bed .list, .picard, .interval_list, .intervals.
     * @param str token to identify as a filename.
     * @return true if the token looks like a filename, or false otherwise.
     */
    public static boolean isIntervalFile(final String str) {
        return isIntervalFile(str, true);
    }

    /**
     * Check if string argument was intended as a file
     * Accepted file extensions are defined in {@link #INTERVAL_FILE_EXTENSIONS}
     * @param str token to identify as a filename.
     * @param checkExists if true throws an exception if the file doesn't exist and has an interval file extension
     * @return true if the token looks like an interval file name, or false otherwise.
     */
    public static boolean isIntervalFile(final String str, final boolean checkExists) {
        Utils.nonNull(str);
        boolean hasIntervalFileExtension = false;
        for ( final String extension : INTERVAL_FILE_EXTENSIONS ) {
            if ( str.toLowerCase().endsWith(extension) ) {
                hasIntervalFileExtension = true;
            }
        }

        if ( hasIntervalFileExtension ) {
            final Path path = IOUtils.getPath(str);
            if ( ! checkExists || Files.exists(path) ) {
                return true;
            } else {
                throw new UserException.CouldNotReadInputFile(path, "The interval file does not exist.");
            }
        }
        else {
            return false;
        }
    }

    /**
     * Returns a map of contig names with their sizes.
     * @param reference The reference for the intervals.
     * @return A map of contig names with their sizes.
     */
    public static Map<String, Integer> getContigSizes(final Path reference) {
        final ReferenceSequenceFile referenceSequenceFile = createReference(reference);
        final List<GenomeLoc> locs = GenomeLocSortedSet.createSetFromSequenceDictionary(referenceSequenceFile.getSequenceDictionary()).toList();
        final Map<String, Integer> lengths = new LinkedHashMap<>();
        for (final GenomeLoc loc: locs) {
            lengths.put(loc.getContig(), loc.size());
        }
        return lengths;
    }

    /**
     * Splits an interval list into multiple files.
     * @param fileHeader The sam file header.
     * @param locs The genome locs to split.
     * @param scatterParts The output interval lists to write to.
     */
    public static void scatterContigIntervals(final SAMFileHeader fileHeader, final List<GenomeLoc> locs, final List<File> scatterParts) {

        // Contract: must divide locs up so that each of scatterParts gets a sublist such that:
        // (a) all locs concerning a particular contig go to the same part
        // (b) locs are not split or combined, and remain in the same order (so scatterParts[0] + ... + scatterParts[n] == locs)

        // Locs are already sorted.

        long totalBases = 0;
        for(final GenomeLoc loc : locs) {
            totalBases += loc.size();
        }

        final long idealBasesPerPart = totalBases / scatterParts.size();
        if(idealBasesPerPart == 0) {
            throw new UserException.BadInput(String.format("Genome region is too short (%d bases) to split into %d parts", totalBases, scatterParts.size()));
        }

        // Find the indices in locs where we switch from one contig to the next.
        final ArrayList<Integer> contigStartLocs = new ArrayList<>();
        String prevContig = null;

        for(int i = 0; i < locs.size(); ++i) {

            final GenomeLoc loc = locs.get(i);
            if(prevContig == null || !loc.getContig().equals(prevContig)) {
                contigStartLocs.add(i);
            }
            prevContig = loc.getContig();

        }

        if(contigStartLocs.size() < scatterParts.size()) {
            throw new UserException.BadInput(String.format("Input genome region has too few contigs (%d) to split into %d parts", contigStartLocs.size(), scatterParts.size()));
        }

        long thisPartBases = 0;
        int partIdx = 0;
        IntervalList outList = new IntervalList(fileHeader);

        for(int i = 0; i < locs.size(); ++i) {

            final GenomeLoc loc = locs.get(i);
            thisPartBases += loc.getStop() - loc.getStart();

            outList.add(toInterval(loc, i));

            boolean partMustStop = false;

            if(partIdx < (scatterParts.size() - 1)) {

                // If there are n contigs and n parts remaining then we must split here,
                // otherwise we will run out of contigs.

                final int nextPart = partIdx + 1;
                final int nextPartMustStartBy = contigStartLocs.get(nextPart + (contigStartLocs.size() - scatterParts.size()));
                if(i + 1 == nextPartMustStartBy) {
                    partMustStop = true;
                }

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
                if((i + 1) < locs.size()) {
                    nextLoc = locs.get(i + 1);
                }

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
    public static List<List<GenomeLoc>> splitIntervalsToSubLists(final List<GenomeLoc> locs, final List<Integer> splits) {
        Utils.nonNull(locs, "locs is null");
        Utils.nonNull(splits, "splits is null");
        int start = 0;
        final List<List<GenomeLoc>> sublists = new ArrayList<>(splits.size());
        for (final Integer stop: splits) {
            final List<GenomeLoc> curList = new ArrayList<>();
            for (int i = start; i < stop; i++) {
                curList.add(locs.get(i));
            }
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
    public static void scatterFixedIntervals(final SAMFileHeader fileHeader, final List<List<GenomeLoc>> splits, final List<File> scatterParts) {
        Utils.nonNull(fileHeader, "fileHeader is null");
        Utils.nonNull(splits, "splits is null");
        Utils.nonNull(scatterParts, "scatterParts is null");
        Utils.containsNoNull(splits, "null split loc");

        if (splits.size() != scatterParts.size()) {
            throw new CommandLineException.BadArgumentValue("splits", String.format("Split points %d does not equal the number of scatter parts %d.", splits.size(), scatterParts.size()));
        }

        int fileIndex = 0;
        int locIndex = 1;
        for (final List<GenomeLoc> split : splits) {
            final IntervalList intervalList = new IntervalList(fileHeader);
            for (final GenomeLoc loc : split) {
                intervalList.add(toInterval(loc, locIndex++));
            }
            intervalList.write(scatterParts.get(fileIndex++));
        }
    }

    /**
     * Splits the genome locs up by size.
     * @param locs Genome locs to split.
     * @param numParts Number of parts to split the locs into.
     * @return The stop points to split the genome locs.
     */
    public static List<List<GenomeLoc>> splitFixedIntervals(final List<GenomeLoc> locs, final int numParts) {
        Utils.nonNull(locs, "locs is null");

        if (locs.size() < numParts) {
            throw new CommandLineException.BadArgumentValue("scatterParts", String.format("Cannot scatter %d locs into %d parts.", locs.size(), numParts));
        }
        final long locsSize = intervalSize(locs);
        final List<Integer> splitPoints = new ArrayList<>();
        addFixedSplit(splitPoints, locs, locsSize, 0, locs.size(), numParts);
        Collections.sort(splitPoints);
        splitPoints.add(locs.size());
        return splitIntervalsToSubLists(locs, splitPoints);
    }

    public static List<List<GenomeLoc>> splitLocusIntervals(final List<GenomeLoc> locs, final int numParts) {
        Utils.nonNull(locs, "locs is null");
        if (numParts < 0) {
            throw new CommandLineException.BadArgumentValue("scatterParts", String.format("Cannot scatter %d locs into %d parts.", locs.size(), numParts));
        }

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

    private static SplitLocusRecursive splitLocusIntervals1(final LinkedList<GenomeLoc> remaining, final long idealSplitSize) {
        final List<GenomeLoc> split = new ArrayList<>();
        long size = 0;

        while ( ! remaining.isEmpty() ) {
            final GenomeLoc head = remaining.pop();
            final long newSize = size + head.size();

            if ( newSize == idealSplitSize ) {
                split.add(head);
                break; // we are done
            } else if ( newSize > idealSplitSize ) {
                final long remainingBp = idealSplitSize - size;
                final long cutPoint = head.getStart() + remainingBp;
                final GenomeLoc[] parts = head.split((int)cutPoint);
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
    public static boolean overlaps(final Locatable left, final Locatable right) {
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

    public static String locatableToString(Locatable interval) {
        Utils.nonNull(interval);
        return String.format("%s:%s-%s", interval.getContig(), interval.getStart(), interval.getEnd());
    }

    private static final class SplitLocusRecursive {
        final List<GenomeLoc> split;
        final LinkedList<GenomeLoc> remaining;

        private SplitLocusRecursive(final List<GenomeLoc> split, final LinkedList<GenomeLoc> remaining) {
            this.split = split;
            this.remaining = remaining;
        }
    }

    public static List<GenomeLoc> flattenSplitIntervals(final List<List<GenomeLoc>> splits) {
        Utils.nonNull(splits, "splits is null");

        final List<GenomeLoc> locs = new ArrayList<>();
        splits.forEach(locs::addAll);
        return locs;
    }

    private static void addFixedSplit(final List<Integer> splitPoints, final List<GenomeLoc> locs, final long locsSize, final int startIndex, final int stopIndex, final int numParts) {
        Utils.nonNull(splitPoints, "splitPoints is null");
        if (numParts < 2) {
            return;
        }
        final int halfParts = (numParts + 1) / 2;
        final Pair<Integer, Long> splitPoint = getFixedSplit(locs, locsSize, startIndex, stopIndex, halfParts, numParts - halfParts);
        final int splitIndex = splitPoint.getLeft();
        final long splitSize = splitPoint.getRight();
        splitPoints.add(splitIndex);
        addFixedSplit(splitPoints, locs, splitSize, startIndex, splitIndex, halfParts);
        addFixedSplit(splitPoints, locs, locsSize - splitSize, splitIndex, stopIndex, numParts - halfParts);
    }

    private static Pair<Integer, Long> getFixedSplit(final List<GenomeLoc> locs, final long locsSize, final int startIndex, final int stopIndex, final int minLocs, final int maxLocs) {
        int splitIndex = startIndex;
        long splitSize = 0;
        for (int i = 0; i < minLocs; i++) {
            splitSize += locs.get(splitIndex).size();
            splitIndex++;
        }
        final long halfSize = locsSize / 2;
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
    private static Interval toInterval(final GenomeLoc loc, final int locIndex) {
        return new Interval(loc.getContig(), loc.getStart(), loc.getStop(), false, "interval_" + locIndex);
    }

    /**
     * merge a list of genome locs that may be overlapping, returning the list of unique genomic locations
     *
     * @param raw the unchecked genome loc list
     * @param rule the merging rule we're using
     *
     * @return the list of merged locations
     */
    public static List<GenomeLoc> mergeIntervalLocations(final List<GenomeLoc> raw, final IntervalMergingRule rule) {
        if (raw.size() <= 1) {
            return Collections.unmodifiableList(raw);
        } else {
            final ArrayList<GenomeLoc> merged = new ArrayList<>();
            final Iterator<GenomeLoc> it = raw.iterator();
            GenomeLoc prev = it.next();
            while (it.hasNext()) {
                final GenomeLoc curr = it.next();
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
        for ( final GenomeLoc loc : locs ) {
            size += loc.size();
        }
        return size;
    }

    /**
     * Returns a list of intervals between the passed int locs. Does not extend UNMAPPED locs.
     * @param parser A genome loc parser for creating the new intervals
     * @param locs Original genome locs
     * @param basePairs Number of base pairs on each side of loc
     * @return The list of intervals between the locs
     */
    public static List<GenomeLoc> getIntervalsWithFlanks(final GenomeLocParser parser, final List<GenomeLoc> locs, final int basePairs) {

        if (locs.isEmpty()) {
            return Collections.emptyList();
        }

        final List<GenomeLoc> expanded = locs.stream()
                .map(loc -> parser.createPaddedGenomeLoc(loc, basePairs))
                .collect(Collectors.toList());

        return sortAndMergeIntervals(parser, expanded, IntervalMergingRule.ALL).toList();
    }

    /**
     * Pads the provided intervals by the specified amount, sorts the resulting intervals, and merges intervals
     * that are adjacent/overlapping after padding.
     *
     * @param intervals intervals to pad
     * @param basePairs number of bases of padding to add to each side of each interval
     * @param dictionary sequence dictionary used to restrict padded intervals to the bounds of their contig
     * @return the provided intervals padded by the specified amount, sorted, with adjacent/overlapping intervals merged
     */
    public static List<SimpleInterval> getIntervalsWithFlanks(final List<SimpleInterval> intervals, final int basePairs, final SAMSequenceDictionary dictionary) {
        final GenomeLocParser parser = new GenomeLocParser(dictionary);
        final List<GenomeLoc> intervalsAsGenomeLocs = genomeLocsFromLocatables(parser, intervals);
        final List<GenomeLoc> paddedGenomeLocs = getIntervalsWithFlanks(parser, intervalsAsGenomeLocs, basePairs);
        return convertGenomeLocsToSimpleIntervals(paddedGenomeLocs);
    }

    /**
     * Accepts a sorted List of intervals, and returns a List of Lists of intervals grouped by contig,
     * one List per contig.
     *
     * @param sortedIntervals sorted List of intervals to group by contig
     * @return A List of Lists of intervals, one List per contig
     */
    public static List<List<SimpleInterval>> groupIntervalsByContig(final List<SimpleInterval> sortedIntervals) {
        final List<List<SimpleInterval>> intervalGroups = new ArrayList<>();
        List<SimpleInterval> currentGroup = new ArrayList<>();
        String currentContig = null;

        for ( final SimpleInterval currentInterval : sortedIntervals ) {
            if ( currentContig != null && ! currentContig.equals(currentInterval.getContig()) ) {
                intervalGroups.add(currentGroup);
                currentGroup = new ArrayList<>();
            }

            currentContig = currentInterval.getContig();
            currentGroup.add(currentInterval);
        }

        if ( ! currentGroup.isEmpty() ) {
            intervalGroups.add(currentGroup);
        }

        return intervalGroups;
    }

    private static ReferenceSequenceFile createReference(final Path fastaPath) {
            return CachingIndexedFastaSequenceFile.checkAndCreate(fastaPath);
    }

    private static LinkedHashMap<String, List<GenomeLoc>> splitByContig(final List<GenomeLoc> sorted) {
        final LinkedHashMap<String, List<GenomeLoc>> splits = new LinkedHashMap<>();
        GenomeLoc last = null;
        List<GenomeLoc> contigLocs = null;
        for (final GenomeLoc loc: sorted) {
            if (GenomeLoc.isUnmapped(loc)) {
                continue;
            }
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
        Utils.nonNull(parser, "the input genome-loc parser cannot be null");

        Utils.nonNull(locatables, "the input locatable collection cannot be null");
        Utils.containsNoNull(locatables, "some element in the locatable input collection is null");

        final List<GenomeLoc> result = locatables.stream().map(parser::createGenomeLoc).collect(Collectors.toList());
        return Collections.unmodifiableList(result);
    }

    /**
     * Builds a list of intervals that cover the whole given sequence.
     */
    public static List<SimpleInterval> getAllIntervalsForReference(final SAMSequenceDictionary sequenceDictionary) {
        return GenomeLocSortedSet.createSetFromSequenceDictionary(sequenceDictionary)
                .stream()
                .map(SimpleInterval::new)
                .collect(Collectors.toList());
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
     * and end positions are not valid
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
        final int lastColonIndex = intervalQueryString.lastIndexOf(SimpleInterval.CONTIG_SEPARATOR);
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
            final int lastDashIndex = intervalQueryString.lastIndexOf(SimpleInterval.START_END_SEPARATOR);
            int startPos;
            int endPos;

            // Try to resolve the suffix as a query against the contig represented by the prefix.
            if (intervalQueryString.endsWith(SimpleInterval.END_OF_CONTIG)) {
                // try to resolve as "prefix:nnn+"
                startPos = SimpleInterval.parsePositionThrowOnFailure(intervalQueryString.substring(lastColonIndex + 1, intervalQueryString.length()-1));
                endPos = prefixSequence.getSequenceLength();
            } else if (lastDashIndex > lastColonIndex) {
                // try to resolve as "prefix:start-end"
                startPos = SimpleInterval.parsePositionThrowOnFailure(intervalQueryString.substring(lastColonIndex + 1, lastDashIndex));
                endPos = SimpleInterval.parsePositionThrowOnFailure(intervalQueryString.substring(lastDashIndex + 1, intervalQueryString.length()));
            } else {
                // finally, try to resolve as "prefix:nnn"
                startPos = SimpleInterval.parsePositionThrowOnFailure(intervalQueryString.substring(lastColonIndex + 1, intervalQueryString.length()));
                endPos = startPos;
            }

            if (SimpleInterval.isValid(prefix, startPos, endPos)) {
                // We've pre-tested to validate the positions, so add this interval. This should never throw.
                resolvedIntervals.add(new SimpleInterval(prefix, startPos, endPos));
            } else {
                // Positions don't appear to be valid, but we don't want to throw if there is any other valid
                // interpretation of the query string
                if (resolvedIntervals.isEmpty()) {
                    // validatePositions throws on validation failure, which is guaranteed if we got here
                    SimpleInterval.validatePositions(prefix, startPos, endPos);
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
     * Create a new interval, bounding start and stop by the start and end of contig
     *
     * This function will return null if start and stop cannot be adjusted in any reasonable way
     * to be on the contig.  For example, if start and stop are both past the end of the contig,
     * there's no way to fix this, and null will be returned.
     *
     * @param contig our contig
     * @param start our start as an arbitrary integer (may be negative, etc)
     * @param stop our stop as an arbitrary integer (may be negative, etc)
     * @param contigLength length of the contig
     * @return a valid interval over contig, or null if a meaningful interval cannot be created
     */
    public static SimpleInterval trimIntervalToContig(final String contig, final int start, final int stop, final int contigLength) {
        Utils.nonNull(contig);
        Utils.validateArg(contigLength >= 1, () -> "contigLength should be at least 1 but was " + contigLength);
        final int boundedStart = Math.max(1, start);
        final int boundedStop = Math.min(contigLength, stop);

        if ( boundedStart > contigLength || boundedStop < 1 ){
            // there's no meaningful way to create this interval, as the start and stop are off the contig
            return null;
        } else {
            return new SimpleInterval(contig, boundedStart, boundedStop);
        }
    }

    /**
     * Determines whether the provided interval is within the bounds of its assigned contig according to the provided dictionary
     *
     * @param interval interval to check
     * @param dictionary dictionary to use to validate contig bounds
     * @return true if the interval's contig exists in the dictionary, and the interval is within its bounds, otherwise false
     */
    public static boolean intervalIsOnDictionaryContig( final SimpleInterval interval, final SAMSequenceDictionary dictionary ) {
        Utils.nonNull(interval);
        Utils.nonNull(dictionary);

        final SAMSequenceRecord contigRecord = dictionary.getSequence(interval.getContig());
        if ( contigRecord == null ) {
            return false;
        }

        return interval.getEnd() <= contigRecord.getSequenceLength();
    }

    //-------------------------------------------------------------------------------------------------
    // Utility code related to the notion of splitting intervals at well-defined boundaries,
    // so that whichever intervals you start from, the resulting shards will line up.

    /**
     *
     * Splits the given input intervals into shards of at most the requested size.
     * The shard boundaries lie at integer multiples of shardSize.
     *
     * chr2:1-200 -> chr2:1-100,chr2:101-200
     *
     * This method will return all intervals in RAM.  If you need a solution that is light on RAM usage, though it
     *  returns an iterator, see {@link org.broadinstitute.hellbender.utils.iterators.ShardedIntervalIterator}
     */
    static public List<SimpleInterval> cutToShards(Iterable<SimpleInterval> intervals, int shardSize) {
        ArrayList<SimpleInterval> ret = new ArrayList<>();
        for (SimpleInterval i : intervals) {
            int beginShard = shardIndex(i.getStart(), shardSize);
            int endShard = shardIndex(i.getEnd(), shardSize);
            if (beginShard==endShard) {
                ret.add(i);
                continue;
            }
            // more than one shard: output begin to end-of-shard, then multiple full shards, then begin-of-shard to end.
            ret.add(new SimpleInterval(i.getContig(), i.getStart(), endOfShard(beginShard, shardSize)));
            for (int shard = beginShard+1; shard<endShard; shard++) {
                ret.add(new SimpleInterval(i.getContig(), beginOfShard(shard, shardSize), endOfShard(shard, shardSize)));
            }
            ret.add(new SimpleInterval(i.getContig(), beginOfShard(endShard, shardSize), i.getEnd()));
        }
        return ret;
    }


    /**
     * number of the shard this offset is in. Shards are numbered starting at zero.
     */
    static public int shardIndex(int oneBasedOffset, int shardSize) {
        return ((oneBasedOffset-1) / shardSize);
    }

    /**
     * first offset in this shard (1-based).
     */
    static public int beginOfShard(int shardIndex, int shardSize) {
        return ((shardIndex) * shardSize) + 1;
    }

    /**
     * last offset in this shard (1-based).
     */
    static public int endOfShard(int shardIndex, int shardSize) {
        return beginOfShard(shardIndex+1, shardSize)-1;
    }


    // (end of shard-related code)

    /**
     * Combine the breakpoints of multiple intervals and return a list of locatables based on the updated breakpoints.
     *
     * Suppose we have two lists of locatables:
     * List 1:
     *
     * <pre>
     *     1  1000  2000
     * </pre>
     *
     * List 2:
     *
     * <pre>
     *     1  500  2500
     *     1  2501  3000
     *     1  4000  5000
     * </pre>
     *
     * The result would be:
     * <pre>
     *     1  500  999
     *     1  1000  2000
     *     1  2001  2500
     *     1  2501  3000
     *     1  4000  5000
     * </pre>
     *
     * Note that start breakpoints will always appear as starts of the resulting intervals.
     *
     * <p>
     * Does not alter the input.
     * </p>
     *
     * Any single list of input locatables containing duplicates or overlapping intervals will throw an exception.
     *
     * Intervals are assumed to include the start and end bases.
     *
     * This method performs all necessary sorting.
     *
     * @param unsortedLocatables1 list of locatables
     * @param unsortedLocatables2 list of locatables
     * @param dictionary Sequence dictionary to base the sort.  The order of contigs/sequences in the dictionary is the order of the sorting here.
     * @return Locatables from the combined breakpoints of unsortedLocatables1 and unsortedLocatables2.  If both inputs are null, return an
     *   empty list.  Please note that returned values are new copies.  If exactly one of the inputs is null, this method
     *   returns a copy of of the non-null input.
     */
    public static <T extends Locatable> List<Locatable> combineAndSortBreakpoints(final List<T> unsortedLocatables1,
                                                                                     final List<T> unsortedLocatables2,
                                                                                     final SAMSequenceDictionary dictionary) {

        Utils.nonNull(dictionary);

        final List<T> locatables1 = sortLocatablesBySequenceDictionary(unsortedLocatables1, dictionary);
        final List<T> locatables2 = sortLocatablesBySequenceDictionary(unsortedLocatables2, dictionary);

        if ((locatables1 == null) && (locatables2 == null)) {
            return Collections.emptyList();
        }

        validateNoOverlappingIntervals(locatables1);
        validateNoOverlappingIntervals(locatables2);

        if (CollectionUtils.isEmpty(locatables1)) {
            return locatables2.stream().map(SimpleInterval::new).collect(Collectors.toList());
        }
        if (CollectionUtils.isEmpty(locatables2)) {
            return locatables1.stream().map(SimpleInterval::new).collect(Collectors.toList());
        }

        final List<Locatable> masterList = new ArrayList<>();
        masterList.addAll(locatables1);
        masterList.addAll(locatables2);

        final Set<String> contigs = masterList.stream()
                .map(Locatable::getContig).collect(Collectors.toSet());

        final Map<String, Set<Pair<Integer, IntervalBreakpointType>>> contigToBreakpoints = contigs.stream()
                .collect(Collectors.toMap(Function.identity(), l -> new HashSet<>()));

        // Populate initialized maps with contigs to the break points.  Also, keep a mapping of contigs
        //  to start breakpoints.
        masterList.forEach(l -> {
            contigToBreakpoints.get(l.getContig()).add(Pair.of(l.getStart(), IntervalBreakpointType.START_BREAKPOINT));
            contigToBreakpoints.get(l.getContig()).add(Pair.of(l.getEnd(), IntervalBreakpointType.END_BREAKPOINT));
        });

        final List<Locatable> result = new ArrayList<>();

        for (final String contig : contigs) {

            // Sort the breakpoints for this contig.  Use the pair structure, since we need to differentiate between a
            //  breakpoint that is a start and a breakpoint that is end, yet have the same position.  This is especially
            //  important for single base intervals.
            final List<Pair<Integer, IntervalBreakpointType>> breakpoints = new ArrayList<>(contigToBreakpoints.get(contig));
            breakpoints.sort((p1, p2) -> {
                final int firstComparison = p1.getLeft().compareTo(p2.getLeft());
                if (firstComparison != 0) {
                    return firstComparison;
                } else {
                    // We want start breakpoints before end breakpoints
                    return p1.getRight().compareTo(p2.getRight());
                }
            });
            int numCurrentStarts = 0;
            int numCurrentEnds = 0;
            for (int i = 0; i < (breakpoints.size() - 1); i++) {
                final int currentBreakpoint = breakpoints.get(i).getLeft();
                final int nextBreakpoint = breakpoints.get(i + 1).getLeft();

                final boolean isCurrentBreakpointStart = breakpoints.get(i).getRight() == IntervalBreakpointType.START_BREAKPOINT;
                final boolean isNextBreakpointStart = breakpoints.get(i + 1).getRight() == IntervalBreakpointType.START_BREAKPOINT;

                // if both breakpoints are starts of intervals, then the result is bp1, bp2-1
                // if both breakpoints are ends of intervals, then the result is bp1+1, bp2
                int start = (!isCurrentBreakpointStart && !isNextBreakpointStart ? currentBreakpoint + 1 : currentBreakpoint);
                int end = (isCurrentBreakpointStart && isNextBreakpointStart ? nextBreakpoint - 1 : nextBreakpoint);

                // If the current breakpoint is an end and the next is a start, then we want to shrink both ends.
                //  Note that this could indicate that we are between intervals for one list of locatables,
                //   but not the other.
                final boolean isBetweenIntervals = !isCurrentBreakpointStart && isNextBreakpointStart;
                if (isBetweenIntervals) {
                    start++;
                    end--;
                }

                // If the next breakpoint is a start and the current is end AND we are not in the middle of an interval,
                //   then we do NOT add the interval at all (we are between intervals in both lists).
                if (isCurrentBreakpointStart) {
                    numCurrentStarts++;
                } else {
                    numCurrentEnds++;
                }

                // We also need to check if start > end, since we could have been between intervals that were one adjacent,
                //   so we don't want to add an interval for that.
                if (((!isBetweenIntervals) || (numCurrentStarts > numCurrentEnds)) &&
                        start <= end){
                    result.add(new SimpleInterval(contig, start, end));
                }
            }
        }
        return sortLocatablesBySequenceDictionary(result, dictionary);
    }

    /**
     *  Throws Bad Input exception if any overlaps are detected within the list of locatables.
     * @param locatables List of locatables to test.  {@code null} will never throw an exception.
     * @param <T> Locatable class
     */
    public static <T extends Locatable> void validateNoOverlappingIntervals(List<T> locatables) {
        // Do not throw an exception for empty or null lists.
        if (CollectionUtils.isEmpty(locatables)) {
            return;
        }

        final HashSet<T> locatablesSet = new HashSet<>(locatables);
        if (locatablesSet.size() != locatables.size()) {
            throw new UserException.BadInput("Duplicate(s) detected in input:  " + locatables.size() + " intervals had " + (locatables.size() - locatablesSet.size()) + " duplicates.");
        }

        final OverlapDetector<T> overlapDetector = OverlapDetector.create(locatables);

        for (final Locatable locatable : locatables) {
            final Set<T> overlaps = overlapDetector.getOverlaps(locatable);
            if (overlaps.size() > 1) {
                throw new UserException.BadInput("Overlap detected in input:  " + locatable + " overlapped " + StringUtils.join(overlaps, ", "));
            }
        }
    }

    /**
     * Sort by the contig then position as specified by the index order in the given sequence dictionary.
     *
     * @param locatables list of locatables.
     * @param dictionary Never {@code null}
     * @param <T> Locatable
     * @return new list that is sorted using the sequence index of the given dictionary.  Returns {@code null} if locatables
     *   is {@code null}.  Instances in the list are not copies of input.
     */
    public static <T extends Locatable> List<T> sortLocatablesBySequenceDictionary(final Collection<T> locatables, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(dictionary);

        final List<T> result = (locatables == null ? null : new ArrayList<>(locatables));

        if (result != null) {
            result.sort(getDictionaryOrderComparator(dictionary));
        }
        return result;
    }

    /**
     * Creates a map of which locatables (keys) overlap the other list of locatables (vals)
     *
     * Input lists will be sorted sorted by the input dictionary.
     *
     *  Within a single input list, segments cannot overlap (or duplicate).  If this occurs, an exception is thrown.
     *
     * No copies of inputs are created.
     *
     * @param keys -- the intervals we wish to query.  Sorted by interval.  No intervals overlap.  Never {@code null}
     * @param vals -- the intervals that we wish to map to the keys.  Sorted by interval.  No intervals overlap.  Never {@code null}
     * @param dictionary -- the SAMSequenceDictionary that the intervals (and sorting) derive from.  Never {@code null}
     * @return a mapping of intervals from keys to the list of overlapping intervals in vals.  All item in keys will
     * have a key.  Never {@code null}
     */
    public static <T extends Locatable, U extends Locatable> Map<T, List<U>> createOverlapMap(final List<T> keys, final List<U> vals,
                                                                                              final SAMSequenceDictionary dictionary) {
        Utils.nonNull(keys);
        Utils.nonNull(vals);
        Utils.nonNull(dictionary);

        validateNoOverlappingIntervals(keys);
        validateNoOverlappingIntervals(vals);

        final List<T> sortedKeys = sortLocatablesBySequenceDictionary(keys, dictionary);

        final OverlapDetector<U> overlapDetector = OverlapDetector.create(vals);

        final Map<T, List<U>> result = new HashMap<>();

        for (final T key: sortedKeys) {

            // Get the overlaps, sort'em, and create a map entry.
            final Set<U> overlaps = overlapDetector.getOverlaps(key);
            final List<U> overlapsAsList = sortLocatablesBySequenceDictionary(overlaps, dictionary);
            result.put(key, overlapsAsList);
        }

        return result;
    }

    /**
     *
     * The order of contigs/sequences in the dictionary is the order of the sorting here.
     *
     * @param dictionary dictionary to use for the sorting.  Intervals with sequences not in this dictionary will cause
     *                   exceptions to be thrown.  Never {@code null}.
     * @return an instance of {@code Comapator<Locatable>} for use in sorting of Locatables.
     */
    public static Comparator<Locatable> getDictionaryOrderComparator(final SAMSequenceDictionary dictionary) {
        Utils.nonNull(dictionary);
        return (o1, o2) -> IntervalUtils.compareLocatables(o1, o2, dictionary);
    }

    /**
     * An enum to classify breakpoints whether the breakpoint is the start or end of a region.
     */
    public enum IntervalBreakpointType {
        START_BREAKPOINT, END_BREAKPOINT
    }
}
