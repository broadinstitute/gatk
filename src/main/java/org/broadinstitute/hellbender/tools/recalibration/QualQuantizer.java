/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2014 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.hellbender.tools.recalibration;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;

import java.io.PrintStream;
import java.util.*;

/**
 * A general algorithm for quantizing quality score distributions to use a specific number of levels
 *
 * Takes a histogram of quality scores and a desired number of levels and produces a
 * map from original quality scores -> quantized quality scores.
 *
 * Note that this data structure is fairly heavy-weight, holding lots of debugging and
 * calculation information.  If you want to use it efficiently at scale with lots of
 * read groups the right way to do this:
 *
 * Map<ReadGroup, List<Byte>> map
 * for each read group rg:
 *   hist = getQualHist(rg)
 *   QualQuantizer qq = new QualQuantizer(hist, nLevels, minInterestingQual)
 *   map.set(rg, qq.getOriginalToQuantizedMap())
 *
 * This map would then be used to look up the appropriate original -> quantized
 * quals for each read as it comes in.
 */
public class QualQuantizer {
    final private static Set<QualInterval> MY_EMPTY_SET = Collections.emptySet();

    private static Logger logger = LogManager.getLogger(QualQuantizer.class);

    /**
     * Inputs to the QualQuantizer
     */
    final int nLevels, minInterestingQual;
    final List<Long> nObservationsPerQual;

    /**
     * Map from original qual (e.g., Q30) to new quantized qual (e.g., Q28).
     *
     * Has the same range as nObservationsPerQual
     */
    final List<Byte> originalToQuantizedMap;

    /** Sorted set of qual intervals.
     *
     * After quantize() this data structure contains only the top-level qual intervals
     */
    final TreeSet<QualInterval> quantizedIntervals;

    /**
     * Protected creator for testng use only
     */
    protected QualQuantizer(final int minInterestingQual) {
        this.nObservationsPerQual = Collections.emptyList();
        this.nLevels = 0;
        this.minInterestingQual = minInterestingQual;
        this.quantizedIntervals = null;
        this.originalToQuantizedMap = null;
    }

    /**
     * Creates a QualQuantizer for the histogram that has nLevels
     *
     * Note this is the only interface to the system.  After creating this object
     * the map can be obtained via getOriginalToQuantizedMap()
     *
     * @param nObservationsPerQual A histogram of counts of bases with quality scores.  Note that
     *  this histogram must start at 0 (i.e., get(0) => count of Q0 bases) and must include counts all the
     *  way up to the largest quality score possible in the reads.  OK if the histogram includes many 0
     *  count bins, as these are quantized for free.
     * @param nLevels the desired number of distinct quality scores to represent the full original range.  Must
     *  be at least 1.
     * @param minInterestingQual All quality scores <= this value are considered uninteresting and are freely
     *  merged together.  For example, if this value is 10, then Q0-Q10 are all considered free to merge, and
     *  quantized into a single value. For ILMN data with lots of Q2 bases this results in a Q2 bin containing
     *  all data with Q0-Q10.
     */
    public QualQuantizer(final List<Long> nObservationsPerQual, final int nLevels, final int minInterestingQual) {
        this.nObservationsPerQual = nObservationsPerQual;
        this.nLevels = nLevels;
        this.minInterestingQual = minInterestingQual;

        // some sanity checking
        if ( Collections.min(nObservationsPerQual) < 0 ) throw new GATKException("Quality score histogram has negative values at: " + Utils.join(", ", nObservationsPerQual));
        if ( nLevels < 0 ) throw new GATKException("nLevels must be >= 0");
        if ( minInterestingQual < 0 ) throw new GATKException("minInterestingQual must be >= 0");

        // actually run the quantizer
        this.quantizedIntervals = quantize();

        // store the map
        this.originalToQuantizedMap = intervalsToMap(quantizedIntervals);
    }

    /**
     * Represents an contiguous interval of quality scores.
     *
     * qStart and qEnd are inclusive, so qStart = qEnd = 2 is the quality score bin of 2
     */
    protected final class QualInterval implements Comparable<QualInterval> {
        final int qStart, qEnd, fixedQual, level;
        final long nObservations, nErrors;
        final Set<QualInterval> subIntervals;

        /** for debugging / visualization.  When was this interval created? */
        int mergeOrder;

        protected QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors, final int level) {
            this(qStart, qEnd, nObservations, nErrors, level, -1, MY_EMPTY_SET);
        }

        protected QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors, final int level, final Set<QualInterval> subIntervals) {
            this(qStart, qEnd, nObservations, nErrors, level, -1, subIntervals);
        }

        protected QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors, final int level, final int fixedQual) {
            this(qStart, qEnd, nObservations, nErrors, level, fixedQual, MY_EMPTY_SET);
        }

        public QualInterval(final int qStart, final int qEnd, final long nObservations, final long nErrors, final int level, final int fixedQual, final Set<QualInterval> subIntervals) {
            this.qStart = qStart;
            this.qEnd = qEnd;
            this.nObservations = nObservations;
            this.nErrors = nErrors;
            this.fixedQual = fixedQual;
            this.level = level;
            this.mergeOrder = 0;
            this.subIntervals = Collections.unmodifiableSet(subIntervals);
        }

        /**
         * @return Human readable name of this interval: e.g., 10-12
         */
        public String getName() {
            return qStart + "-" + qEnd;
        }

        @Override
        public String toString() {
            return "QQ:" + getName();
        }

        /**
         * @return the error rate (in real space) of this interval, or 0 if there are no observations
         */
        public double getErrorRate() {
            if ( hasFixedQual() )
                return QualityUtils.qualToErrorProb((byte) fixedQual);
            else if ( nObservations == 0 )
                return 0.0;
            else
                return (nErrors+1) / (1.0 * (nObservations+1));
        }

        /**
         * @return the QUAL of the error rate of this interval, or the fixed qual if this interval was created with a fixed qual.
         */
        public byte getQual() {
            if ( ! hasFixedQual() )
                return QualityUtils.errorProbToQual(getErrorRate());
            else
                return (byte)fixedQual;
        }

        /**
         * @return true if this bin is using a fixed qual
         */
        public boolean hasFixedQual() {
            return fixedQual != -1;
        }

        @Override
        public int compareTo(final QualInterval qualInterval) {
            return Integer.valueOf(this.qStart).compareTo(qualInterval.qStart);
        }

        /**
         * Create a interval representing the merge of this interval and toMerge
         *
         * Errors and observations are combined
         * Subintervals updated in order of left to right (determined by qStart)
         * Level is 1 + highest level of this and toMerge
         * Order must be updated elsewhere
         *
         * @param toMerge
         * @return newly created merged QualInterval
         */
        public QualInterval merge(final QualInterval toMerge) {
            final QualInterval left = this.compareTo(toMerge) < 0 ? this : toMerge;
            final QualInterval right = this.compareTo(toMerge) < 0 ? toMerge : this;

            if ( left.qEnd + 1 != right.qStart )
                throw new GATKException("Attempting to merge non-contiguous intervals: left = " + left + " right = " + right);

            final long nCombinedObs = left.nObservations + right.nObservations;
            final long nCombinedErr = left.nErrors + right.nErrors;

            final int level = Math.max(left.level, right.level) + 1;
            final Set<QualInterval> subIntervals = new HashSet<QualInterval>(Arrays.asList(left, right));
            QualInterval merged = new QualInterval(left.qStart, right.qEnd, nCombinedObs, nCombinedErr, level, subIntervals);

            return merged;
        }

        public double getPenalty() {
            return calcPenalty(getErrorRate());
        }


        /**
         * Calculate the penalty of this interval, given the overall error rate for the interval
         *
         * If the globalErrorRate is e, this value is:
         *
         * sum_i |log10(e_i) - log10(e)| * nObservations_i
         *
         * each the index i applies to all leaves of the tree accessible from this interval
         * (found recursively from subIntervals as necessary)
         *
         * @param globalErrorRate overall error rate in real space against which we calculate the penalty
         * @return the cost of approximating the bins in this interval with the globalErrorRate
         */
        private double calcPenalty(final double globalErrorRate) {
            if ( globalErrorRate == 0.0 ) // there were no observations, so there's no penalty
                return 0.0;

            if ( subIntervals.isEmpty() ) {
                // this is leave node
                if ( this.qEnd <= minInterestingQual )
                    // It's free to merge up quality scores below the smallest interesting one
                    return 0;
                else {
                    return (Math.abs(Math.log10(getErrorRate()) - Math.log10(globalErrorRate))) * nObservations;
                }
            } else {
                double sum = 0;
                for ( final QualInterval interval : subIntervals )
                    sum += interval.calcPenalty(globalErrorRate);
                return sum;
            }
        }
    }

    /**
     * Main method for computing the quantization intervals.
     *
     * Invoked in the constructor after all input variables are initialized.  Walks
     * over the inputs and builds the min. penalty forest of intervals with exactly nLevel
     * root nodes.  Finds this min. penalty forest via greedy search, so is not guarenteed
     * to find the optimal combination.
     *
     * TODO: develop a smarter algorithm
     *
     * @return the forest of intervals with size == nLevels
     */
    private TreeSet<QualInterval> quantize() {
        // create intervals for each qual individually
        final TreeSet<QualInterval> intervals = new TreeSet<QualInterval>();
        for ( int qStart = 0; qStart < getNQualsInHistogram(); qStart++ ) {
            final long nObs = nObservationsPerQual.get(qStart);
            final double errorRate = QualityUtils.qualToErrorProb((byte)qStart);
            final double nErrors = nObs * errorRate;
            final QualInterval qi = new QualInterval(qStart, qStart, nObs, (int) Math.floor(nErrors), 0, (byte)qStart);
            intervals.add(qi);
        }

        // greedy algorithm:
        // while ( n intervals >= nLevels ):
        //   find intervals to merge with least penalty
        //   merge it
        while ( intervals.size() > nLevels ) {
            mergeLowestPenaltyIntervals(intervals);
        }

        return intervals;
    }

    /**
     * Helper function that finds and merges together the lowest penalty pair of intervals
     * @param intervals
     */
    private void mergeLowestPenaltyIntervals(final TreeSet<QualInterval> intervals) {
        // setup the iterators
        final Iterator<QualInterval> it1 = intervals.iterator();
        final Iterator<QualInterval> it1p = intervals.iterator();
        it1p.next(); // skip one

        // walk over the pairs of left and right, keeping track of the pair with the lowest merge penalty
        QualInterval minMerge = null;
        if ( logger.isDebugEnabled() ) logger.debug("mergeLowestPenaltyIntervals: " + intervals.size());
        int lastMergeOrder = 0;
        while ( it1p.hasNext() ) {
            final QualInterval left = it1.next();
            final QualInterval right = it1p.next();
            final QualInterval merged = left.merge(right);
            lastMergeOrder = Math.max(Math.max(lastMergeOrder, left.mergeOrder), right.mergeOrder);
            if ( minMerge == null || (merged.getPenalty() < minMerge.getPenalty() ) ) {
                if ( logger.isDebugEnabled() ) logger.debug("  Updating merge " + minMerge);
                minMerge = merged;
            }
        }

        // now actually go ahead and merge the minMerge pair
        if ( logger.isDebugEnabled() ) logger.debug("  => final min merge " + minMerge);
        intervals.removeAll(minMerge.subIntervals);
        intervals.add(minMerge);
        minMerge.mergeOrder = lastMergeOrder + 1;
        if ( logger.isDebugEnabled() ) logger.debug("updated intervals: " + intervals);
    }

    /**
     * Given a final forest of intervals constructs a list mapping
     * list.get(i) => quantized qual to use for original quality score i
     *
     * This function should be called only once to initialize the corresponding
     * cached value in this object, as the calculation is a bit costly.
     *
     * @param intervals
     * @return
     */
    private List<Byte> intervalsToMap(final TreeSet<QualInterval> intervals) {
        final List<Byte> map = new ArrayList<Byte>(getNQualsInHistogram());
        map.addAll(Collections.nCopies(getNQualsInHistogram(), Byte.MIN_VALUE));
        for ( final QualInterval interval : intervals ) {
            for ( int q = interval.qStart; q <= interval.qEnd; q++ ) {
                map.set(q, interval.getQual());
            }
        }

        if ( Collections.min(map) == Byte.MIN_VALUE )
            throw new GATKException("quantized quality score map contains an un-initialized value");

        return map;
    }

    private final int getNQualsInHistogram() {
        return nObservationsPerQual.size();
    }

    /**
     * Write out a GATKReport to visualize the QualQuantization process of this data
     * @param out
     */
    public void writeReport(PrintStream out) {
        final GATKReport report = new GATKReport();

        addQualHistogramToReport(report);
        addIntervalsToReport(report);

        report.print(out);
    }

    private final void addQualHistogramToReport(final GATKReport report) {
        report.addTable("QualHistogram", "Quality score histogram provided to report", 2);
        GATKReportTable table = report.getTable("QualHistogram");

        table.addColumn("qual");
        table.addColumn("count");

        for ( int q = 0; q < nObservationsPerQual.size(); q++ ) {
            table.set(q, "qual", q);
            table.set(q, "count", nObservationsPerQual.get(q));
        }
    }


    private final void addIntervalsToReport(final GATKReport report) {
        report.addTable("QualQuantizerIntervals", "Table of QualQuantizer quantization intervals", 10);
        GATKReportTable table = report.getTable("QualQuantizerIntervals");

        table.addColumn("name");
        table.addColumn("qStart");
        table.addColumn("qEnd");
        table.addColumn("level");
        table.addColumn("merge.order");
        table.addColumn("nErrors");
        table.addColumn("nObservations");
        table.addColumn("qual");
        table.addColumn("penalty");
        table.addColumn("root.node");
        //table.addColumn("subintervals", "NA");

        for ( QualInterval interval : quantizedIntervals )
            addIntervalToReport(table, interval, true);
    }

    private final void addIntervalToReport(final GATKReportTable table, final QualInterval interval, final boolean atRootP) {
        final String name = interval.getName();
        table.set(name, "name", name);
        table.set(name, "qStart", interval.qStart);
        table.set(name, "qEnd", interval.qEnd);
        table.set(name, "level", interval.level);
        table.set(name, "merge.order", interval.mergeOrder);
        table.set(name, "nErrors", interval.nErrors);
        table.set(name, "nObservations", interval.nObservations);
        table.set(name, "qual", interval.getQual());
        table.set(name, "penalty", String.format("%.1f", interval.getPenalty()));
        table.set(name, "root.node", atRootP);

        for ( final QualInterval sub : interval.subIntervals )
            addIntervalToReport(table, sub, false);
    }

    public List<Byte> getOriginalToQuantizedMap() {
        return originalToQuantizedMap;
    }
}
