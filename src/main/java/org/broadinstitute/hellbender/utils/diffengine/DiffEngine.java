/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.diffengine;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.classloader.PluginManager;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.report.GATKReport;
import org.broadinstitute.gatk.utils.report.GATKReportTable;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/4/11
 * Time: 12:51 PM
 * A generic engine for comparing tree-structured objects
 *
 */
public class DiffEngine {
    final protected static Logger logger = Logger.getLogger(DiffEngine.class);

    private final Map<String, DiffableReader> readers = new HashMap<String, DiffableReader>();

    public DiffEngine() {
        loadDiffableReaders();
    }

    // --------------------------------------------------------------------------------
    //
    // difference calculation
    //
    // --------------------------------------------------------------------------------

    public List<Difference> diff(DiffElement master, DiffElement test) {
        DiffValue masterValue = master.getValue();
        DiffValue testValue = test.getValue();

        if ( masterValue.isCompound() && masterValue.isCompound() ) {
            return diff(master.getValueAsNode(), test.getValueAsNode());
        } else if ( masterValue.isAtomic() && testValue.isAtomic() ) {
            return diff(masterValue, testValue);
        } else {
            // structural difference in types.  one is node, other is leaf
            return Arrays.asList(new Difference(master, test));
        }
    }

    public List<Difference> diff(DiffNode master, DiffNode test) {
        Set<String> allNames = new HashSet<String>(master.getElementNames());
        allNames.addAll(test.getElementNames());
        List<Difference> diffs = new ArrayList<Difference>();

        for ( String name : allNames ) {
            DiffElement masterElt = master.getElement(name);
            DiffElement testElt = test.getElement(name);
            if ( masterElt == null && testElt == null ) {
                throw new ReviewedGATKException("BUG: unexpectedly got two null elements for field: " + name);
            } else if ( masterElt == null || testElt == null ) { // if either is null, we are missing a value
                // todo -- should one of these be a special MISSING item?
                diffs.add(new Difference(masterElt, testElt));
            } else {
                diffs.addAll(diff(masterElt, testElt));
            }
        }

        return diffs;
    }

    public List<Difference> diff(DiffValue master, DiffValue test) {
        if ( master.getValue().equals(test.getValue()) ) {
            return Collections.emptyList();
        } else {
            return Arrays.asList(new Difference(master.getBinding(), test.getBinding()));
        }
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
    public void reportSummarizedDifferences(List<Difference> diffs, SummaryReportParams params ) {
        printSummaryReport(summarizedDifferencesOfPaths(diffs, params.doPairwise, params.maxRawDiffsToSummarize), params );
    }

    final protected static String[] diffNameToPath(String diffName) {
        return diffName.split("\\.");
    }

    protected List<Difference> summarizedDifferencesOfPathsFromString(List<String> singletonDiffs) {
        List<Difference> diffs = new ArrayList<Difference>();

        for ( String diff : singletonDiffs ) {
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
    private Map<String, Difference> initialPairwiseSummaries(final List<? extends Difference> singletonDiffs,
                                                             final int maxRawDiffsToSummarize) {
        Map<String, Difference> summaries = new HashMap<String, Difference>();

        // create the initial set of differences
        for ( int i = 0; i < singletonDiffs.size(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                Difference diffPath1 = singletonDiffs.get(i);
                Difference diffPath2 = singletonDiffs.get(j);
                if ( diffPath1.length() == diffPath2.length() ) {
                    int lcp = longestCommonPostfix(diffPath1.getParts(), diffPath2.getParts());
                    String path = diffPath2.getPath();
                    if ( lcp != 0 && lcp != diffPath1.length() )
                        path = summarizedPath(diffPath2.getParts(), lcp);
                    Difference sumDiff = new Difference(path, diffPath2.getMaster(), diffPath2.getTest());
                    sumDiff.setCount(0);
                    addSummaryIfMissing(summaries, sumDiff);

                    if ( maxRawDiffsToSummarize != -1 && summaries.size() > maxRawDiffsToSummarize)
                        return summaries;
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
    private Map<String, Difference> initialLeafSummaries(final List<? extends Difference> singletonDiffs,
                                                         final int maxRawDiffsToSummarize) {
        Map<String, Difference> summaries = new HashMap<String, Difference>();

        // create the initial set of differences
        for ( final Difference d : singletonDiffs ) {
            final String path = summarizedPath(d.getParts(), 1);
            Difference sumDiff = new Difference(path, d.getMaster(), d.getTest());
            sumDiff.setCount(0);
            addSummaryIfMissing(summaries, sumDiff);

            if ( maxRawDiffsToSummarize != -1 && summaries.size() > maxRawDiffsToSummarize)
                return summaries;
        }

        return summaries;
    }

    protected List<Difference> summarizedDifferencesOfPaths(final List<? extends Difference> singletonDiffs,
                                                            final boolean doPairwise,
                                                            final int maxRawDiffsToSummarize) {
        final Map<String, Difference> summaries = doPairwise
                ? initialPairwiseSummaries(singletonDiffs, maxRawDiffsToSummarize)
                : initialLeafSummaries(singletonDiffs, maxRawDiffsToSummarize);

        // count differences
        for ( Difference diffPath : singletonDiffs ) {
            for ( Difference sumDiff : summaries.values() ) {
                if ( sumDiff.matches(diffPath.getParts()) )
                    sumDiff.incCount();
            }
        }

        List<Difference> sortedSummaries = new ArrayList<Difference>(summaries.values());
        Collections.sort(sortedSummaries);
        return sortedSummaries;
    }

    protected void addSummaryIfMissing(Map<String, Difference> summaries, Difference diff) {
        if ( ! summaries.containsKey(diff.getPath()) ) {
            summaries.put(diff.getPath(), diff);
        }
    }

    protected void printSummaryReport(List<Difference> sortedSummaries, SummaryReportParams params ) {
        List<Difference> toShow = new ArrayList<Difference>();
        int count = 0, count1 = 0;
        for ( Difference diff : sortedSummaries ) {
            if ( diff.getCount() < params.minSumDiffToShow )
                // in order, so break as soon as the count is too low
                break;

            if ( params.maxItemsToDisplay != 0 && count++ > params.maxItemsToDisplay )
                break;

            if ( diff.getCount() == 1 ) {
                count1++;
                if ( params.maxCountOneItems != 0 && count1 > params.maxCountOneItems )
                    break;
            }

            toShow.add(diff);
        }

        // if we want it in descending order, reverse the list
        if ( ! params.descending ) {
            Collections.reverse(toShow);
        }

        // now that we have a specific list of values we want to show, display them
        GATKReport report = new GATKReport();
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
        GATKReport output = new GATKReport(table);
        output.print(params.out);
    }

    protected static int longestCommonPostfix(String[] diffPath1, String[] diffPath2) {
        int i = 0;
        for ( ; i < diffPath1.length; i++ ) {
            int j = diffPath1.length - i - 1;
            if ( ! diffPath1[j].equals(diffPath2[j]) )
                break;
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
    protected static String summarizedPath(String[] parts, int commonPostfixLength) {
        int stop = parts.length - commonPostfixLength;
        if ( stop > 0 ) parts = parts.clone();
        for ( int i = 0; i < stop; i++ ) {
            parts[i] = "*";
        }
        return Utils.join(".", parts);
    }

    // --------------------------------------------------------------------------------
    //
    // plugin manager
    //
    // --------------------------------------------------------------------------------

    public void loadDiffableReaders() {
        List<Class<? extends DiffableReader>> drClasses = new PluginManager<DiffableReader>( DiffableReader.class ).getPlugins();

        logger.info("Loading diffable modules:");
        for (Class<? extends DiffableReader> drClass : drClasses ) {
            logger.info("\t" + drClass.getSimpleName());

            try {
                DiffableReader dr = drClass.newInstance();
                readers.put(dr.getName(), dr);
            } catch (InstantiationException e) {
                throw new ReviewedGATKException("Unable to instantiate module '" + drClass.getSimpleName() + "'");
            } catch (IllegalAccessException e) {
                throw new ReviewedGATKException("Illegal access error when trying to instantiate '" + drClass.getSimpleName() + "'");
            }
        }
    }

    protected Map<String, DiffableReader> getReaders() {
        return readers;
    }

    protected DiffableReader getReader(String name) {
        return readers.get(name);
    }

    /**
     * Returns a reader appropriate for this file, or null if no such reader exists
     * @param file
     * @return
     */
    public DiffableReader findReaderForFile(File file) {
        for ( DiffableReader reader : readers.values() )
            if (reader.canRead(file) )
                return reader;

        return null;
    }

    /**
     * Returns true if reader appropriate for this file, or false if no such reader exists
     * @param file
     * @return
     */
    public boolean canRead(File file) {
        return findReaderForFile(file) != null;
    }


    public DiffElement createDiffableFromFile(File file) {
        return createDiffableFromFile(file, -1);
    }

    public DiffElement createDiffableFromFile(File file, int maxElementsToRead) {
        DiffableReader reader = findReaderForFile(file);
        if ( reader == null )
            throw new UserException("Unsupported file type: " + file);
        else
            return reader.readFromFile(file, maxElementsToRead);
    }

    public static boolean simpleDiffFiles(File masterFile, File testFile, int maxElementsToRead, DiffEngine.SummaryReportParams params) {
        DiffEngine diffEngine = new DiffEngine();

        if ( diffEngine.canRead(masterFile) && diffEngine.canRead(testFile) ) {
            DiffElement master = diffEngine.createDiffableFromFile(masterFile, maxElementsToRead);
            DiffElement test = diffEngine.createDiffableFromFile(testFile, maxElementsToRead);
            List<Difference> diffs = diffEngine.diff(master, test);
            diffEngine.reportSummarizedDifferences(diffs, params);
            return true;
        } else {
            return false;
        }
    }

    public static class SummaryReportParams {
        final PrintStream out;
        final int maxItemsToDisplay;
        final int maxCountOneItems;
        final int minSumDiffToShow;
        final int maxRawDiffsToSummarize;
        final boolean doPairwise;
        boolean descending = true;

        public SummaryReportParams(PrintStream out,
                                   int maxItemsToDisplay,
                                   int maxCountOneItems,
                                   int minSumDiffToShow,
                                   int maxRawDiffsToSummarize,
                                   final boolean doPairwise) {
            this.out = out;
            this.maxItemsToDisplay = maxItemsToDisplay;
            this.maxCountOneItems = maxCountOneItems;
            this.minSumDiffToShow = minSumDiffToShow;
            this.maxRawDiffsToSummarize = maxRawDiffsToSummarize;
            this.doPairwise = doPairwise;
        }

        public void setDescending(boolean descending) {
            this.descending = descending;
        }
    }
}
