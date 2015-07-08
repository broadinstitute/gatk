package org.broadinstitute.hellbender.utils.diffengine;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;
import org.broadinstitute.hellbender.utils.test.BaseTest;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * A generic engine for comparing tree-structured objects.
 */
public final class DiffEngine {

    private static final int MAX_RECORDS_TO_READ = 1000000;
    private static final int MAX_RAW_DIFFS_TO_SUMMARIZE = -1;

    //Note: diff engines are all equal but static collections are problematic
    //so we create a new diff engine when we need one.

    private final Map<String, DiffableReader> readers = new HashMap<>();

    public DiffEngine() {
        //Make readers for all known diffable formats
        final DiffableReader bamdr = new BAMDiffableReader();
        readers.put(bamdr.getName(), bamdr);

        final DiffableReader vcfdr = new VCFDiffableReader();
        readers.put(vcfdr.getName(), vcfdr);
    }

    // --------------------------------------------------------------------------------
    //
    // Summarizing differences
    //
    // --------------------------------------------------------------------------------

    /**
     * Emits a summary of the diffs to out.  Suppose you have the following three differences:
     *
     *   A.X.Z:1!=2
     *   A.Y.Z:3!=4
     *   B.X.Z:5!=6
     *
     * The above is the itemized list of the differences.  The summary looks for common differences
     * in the name hierarchy, counts those shared elements, and emits the differences that occur
     * in order of decreasing counts.
     *
     * So, in the above example, what are the shared elements?
     *
     * A.X.Z and B.X.Z share X.Z, so there's a *.X.Z with count 2
     * A.X.Z, A.Y.Z, and B.X.Z all share *.*.Z, with count 3
     * Each of A.X.Z, A.Y.Z, and B.X.Z are individually unique, with count 1
     *
     * So we would emit the following summary:
     *
     * *.*.Z: 3
     * *.X.Z: 2
     * A.X.Z: 1 [specific difference: 1!=2]
     * A.Y.Z: 1 [specific difference: 3!=4]
     * B.X.Z: 1 [specific difference: 5!=6]
     *
     * The algorithm to accomplish this calculation is relatively simple. Start with all of the
     * concrete differences.  For each pair of differences A1.A2....AN and B1.B2....BN:
     *
     * find the longest common subsequence Si.Si+1...SN where Ai = Bi = Si
     * If i == 0, then there's no shared substructure
     * If i > 0, then generate the summarized value X = *.*...Si.Si+1...SN
     * if X is a known summary, increment it's count, otherwise set its count to 1
     *
     * Not that only pairs of the same length are considered as potentially equivalent
     *
     * @param params determines how we display the items
     * @param diffs the list of differences to summarize
     */
    static void reportSummarizedDifferences(final List<Difference> diffs, final SummaryReportParams params ) {
        printSummaryReport(summarizedDifferencesOfPaths(diffs, params.doPairwise, params.maxRawDiffsToSummarize), params );
    }

    static String[] diffNameToPath(final String diffName) {
        return diffName.split("\\.");
    }

    static List<Difference> summarizedDifferencesOfPathsFromString(final List<String> singletonDiffs) {
        final List<Difference> diffs = new ArrayList<>();

        for ( final String diff : singletonDiffs ) {
            diffs.add(new Difference(diff));
        }

        return summarizedDifferencesOfPaths(diffs, true, -1);
    }

    /**
     * Computes a minimum set of potential differences between all singleton differences
     * in singletonDiffs.  Employs an expensive pairwise O(n^2) algorithm.
     *
     * @param singletonDiffs
     * @param maxRawDiffsToSummarize
     * @return
     */
    private static Map<String, Difference> initialPairwiseSummaries(final List<? extends Difference> singletonDiffs,
                                                             final int maxRawDiffsToSummarize) {
        final Map<String, Difference> summaries = new HashMap<>();

        // create the initial set of differences
        for ( int i = 0; i < singletonDiffs.size(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                final Difference diffPath1 = singletonDiffs.get(i);
                final Difference diffPath2 = singletonDiffs.get(j);
                if ( diffPath1.length() == diffPath2.length() ) {
                    final int lcp = longestCommonPostfix(diffPath1.getParts(), diffPath2.getParts());
                    String path = diffPath2.getPath();
                    if ( lcp != 0 && lcp != diffPath1.length() ) {
                        path = summarizedPath(diffPath2.getParts(), lcp);
                    }
                    final Difference sumDiff = new Difference(path, diffPath2.getMaster(), diffPath2.getTest());
                    sumDiff.resetCount();;
                    addSummaryIfMissing(summaries, sumDiff);

                    if ( maxRawDiffsToSummarize != -1 && summaries.size() > maxRawDiffsToSummarize) {
                        return summaries;
                    }
                }
            }
        }

        return summaries;
    }

    /**
     * Computes the possible leaf differences among the singleton diffs.
     *
     * The leaf differences are all of the form *.*...*.X where all internal
     * differences are wildcards and the only summarized difference considered
     * interesting to compute is
     *
     * @param singletonDiffs
     * @param maxRawDiffsToSummarize
     * @return
     */
    private static Map<String, Difference> initialLeafSummaries(final List<? extends Difference> singletonDiffs,
                                                         final int maxRawDiffsToSummarize) {
        final Map<String, Difference> summaries = new HashMap<>();

        // create the initial set of differences
        for ( final Difference d : singletonDiffs ) {
            final String path = summarizedPath(d.getParts(), 1);
            final Difference sumDiff = new Difference(path, d.getMaster(), d.getTest());
            sumDiff.resetCount();
            addSummaryIfMissing(summaries, sumDiff);

            if ( maxRawDiffsToSummarize != -1 && summaries.size() > maxRawDiffsToSummarize) {
                return summaries;
            }
        }

        return summaries;
    }

    private static List<Difference> summarizedDifferencesOfPaths(final List<? extends Difference> singletonDiffs,
                                                            final boolean doPairwise,
                                                            final int maxRawDiffsToSummarize) {
        final Map<String, Difference> summaries = doPairwise
                ? initialPairwiseSummaries(singletonDiffs, maxRawDiffsToSummarize)
                : initialLeafSummaries(singletonDiffs, maxRawDiffsToSummarize);

        // count differences
        for ( final Difference diffPath : singletonDiffs ) {
            for ( final Difference sumDiff : summaries.values() ) {
                if ( sumDiff.matches(diffPath.getParts()) ) {
                    sumDiff.incCount();
                }
            }
        }

        final List<Difference> sortedSummaries = new ArrayList<>(summaries.values());
        Collections.sort(sortedSummaries);
        return sortedSummaries;
    }

    private static void addSummaryIfMissing(final Map<String, Difference> summaries, final Difference diff) {
        if ( ! summaries.containsKey(diff.getPath()) ) {
            summaries.put(diff.getPath(), diff);
        }
    }

    private static void printSummaryReport(final List<Difference> sortedSummaries, final SummaryReportParams params ) {
        final List<Difference> toShow = new ArrayList<>();
        int count = 0, count1 = 0;
        for ( final Difference diff : sortedSummaries ) {
            if ( diff.getCount() < params.minSumDiffToShow ) {
                // in order, so break as soon as the count is too low
                break;
            }

            if ( params.maxItemsToDisplay != 0 && count++ > params.maxItemsToDisplay ) {
                break;
            }

            if ( diff.getCount() == 1 ) {
                count1++;
                if ( params.maxCountOneItems != 0 && count1 > params.maxCountOneItems ) {
                    break;
                }
            }

            toShow.add(diff);
        }

        // if we want it in descending order, reverse the list
        if ( ! params.descending ) {
            Collections.reverse(toShow);
        }

        // now that we have a specific list of values we want to show, display them
        final GATKReport report = new GATKReport();
        final String tableName = "differences";
        report.addTable(tableName, "Summarized differences between the master and test files. See http://www.broadinstitute.org/gatk/guide/article?id=1299 for more information", 3);
        final GATKReportTable table = report.getTable(tableName);
        table.addColumn("Difference");
        table.addColumn("NumberOfOccurrences");
        table.addColumn("ExampleDifference");
        for ( final Difference diff : toShow ) {
            final String key = diff.getPath();
            table.addRowID(key, true);
            table.set(key, "NumberOfOccurrences", diff.getCount());
            table.set(key, "ExampleDifference", diff.valueDiffString());
        }
        final GATKReport output = new GATKReport(table);
        output.print(params.out);
    }

    static int longestCommonPostfix(final String[] diffPath1, final String[] diffPath2) {
        int i = 0;
        for ( ; i < diffPath1.length; i++ ) {
            final int j = diffPath1.length - i - 1;
            if ( ! diffPath1[j].equals(diffPath2[j]) ) {
                break;
            }
        }
        return i;
    }

    /**
     * parts is [A B C D]
     * commonPostfixLength: how many parts are shared at the end, suppose its 2
     * We want to create a string *.*.C.D
     *
     * @param parts the separated path values [above without .]
     * @param commonPostfixLength
     * @return
     */
    static String summarizedPath(final String[] partsArg, final int commonPostfixLength) {
        final int stop = partsArg.length - commonPostfixLength;
        if (stop <= 0) {
            return String.join(".", partsArg);
        } else {
            final String[] parts = Arrays.copyOf(partsArg, partsArg.length);
            for ( int i = 0; i < stop; i++ ) {
                parts[i] = "*";
            }
            return String.join(".", parts);
        }
    }


    /**
     * Returns an unmodifiable map of reader names to readers.
     */
    Map<String, DiffableReader> getReaders() {
        return Collections.unmodifiableMap(readers);
    }

    DiffableReader getReader(final String name) {
        return readers.get(name);
    }

    /**
     * Returns a reader appropriate for this file, or null if no such reader exists
     * @param file
     * @return
     */
    DiffableReader findReaderForFile(final File file) {
        for ( final DiffableReader reader : readers.values() ) {
            if (reader.canRead(file)) {
                return reader;
            }
        }

        return null;
    }

    /**
     * Returns true if reader appropriate for this file, or false if no such reader exists
     * @param file
     * @return
     */
    public boolean canRead(final File file) {
        return findReaderForFile(file) != null;
    }

    DiffElement createDiffableFromFile(final File file, final int maxElementsToRead) throws IOException {
        final DiffableReader reader = findReaderForFile(file);
        if ( reader == null )
            throw new UserException("Unsupported file type: " + file);
        else
            return reader.readFromFile(file, maxElementsToRead);
    }

    private static boolean simpleDiffFiles(final File masterFile, final File testFile, final int maxElementsToRead, final DiffEngine.SummaryReportParams params) throws IOException {
        final DiffEngine diffEngine = new DiffEngine();

        if ( diffEngine.canRead(masterFile) && diffEngine.canRead(testFile) ) {
            final DiffElement master = diffEngine.createDiffableFromFile(masterFile, maxElementsToRead);
            final DiffElement test = diffEngine.createDiffableFromFile(testFile, maxElementsToRead);
            final List<Difference> diffs = master.diffElements(test);
            DiffEngine.reportSummarizedDifferences(diffs, params);
            return true;
        } else {
            return false;
        }
    }

    public String diffFiles(final File masterFile, final File testFile) throws IOException {
        String diffEngineOutput = "";
        final ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try(final PrintStream ps = new PrintStream(baos)) {
            final DiffEngine.SummaryReportParams params = new DiffEngine.SummaryReportParams(ps, 20, 10, 0, MAX_RAW_DIFFS_TO_SUMMARIZE, false);
            final boolean success = DiffEngine.simpleDiffFiles(masterFile, testFile, MAX_RECORDS_TO_READ, params);
            if (success) {
                diffEngineOutput = baos.toString();
                BaseTest.log(diffEngineOutput);
                System.out.printf("Note that the above list is not comprehensive.  At most 20 lines of output, and 10 specific differences will be listed.");
            } else {
                return null;
            }
        }
        return diffEngineOutput;
    }

    static final class SummaryReportParams {
        final PrintStream out;
        final int maxItemsToDisplay;
        final int maxCountOneItems;
        final int minSumDiffToShow;
        final int maxRawDiffsToSummarize;
        final boolean doPairwise;
        final boolean descending = true;

        public SummaryReportParams(final PrintStream out,
                                   final int maxItemsToDisplay,
                                   final int maxCountOneItems,
                                   final int minSumDiffToShow,
                                   final int maxRawDiffsToSummarize,
                                   final boolean doPairwise) {
            this.out = out;
            this.maxItemsToDisplay = maxItemsToDisplay;
            this.maxCountOneItems = maxCountOneItems;
            this.minSumDiffToShow = minSumDiffToShow;
            this.maxRawDiffsToSummarize = maxRawDiffsToSummarize;
            this.doPairwise = doPairwise;
        }
    }
}
