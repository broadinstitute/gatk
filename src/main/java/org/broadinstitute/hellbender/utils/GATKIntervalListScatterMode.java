/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import picard.util.IntervalList.IntervalListScatterer;
import picard.util.IntervalList.IntervalListScattererByIntervalCount;
import picard.util.IntervalList.IntervalListScattererWithSubdivision;
import picard.util.IntervalList.IntervalListScattererWithoutSubdivision;
import picard.util.IntervalList.IntervalListScattererWithoutSubdivisionWithOverflow;

import java.util.function.Supplier;

/**
 * An enum to control the creation of the various IntervalListScatter objects
 */
public enum GATKIntervalListScatterMode implements CommandLineParser.ClpEnum {
    /**
     * A simple scatter approach in which all output intervals have size equal to the total base count of the source list divide by the
     * scatter count (except for possible variance in the final interval list).
     */
    INTERVAL_SUBDIVISION(IntervalListScattererWithSubdivision::new, "Scatter the interval list into similarly sized interval lists " +
            "(by base count), breaking up intervals as needed.") {

    },

    /**
     * A scatter approach which attempts to break the intervals into approximately equally weighted pieces where the weights are given as pat
     * of WeightedIntervals.
     */
    INTERVAL_SUBDIVISION_BY_WEIGHT(IntervalListScattererWithSubdivisionByWeight::new, "Scatter the interval list into similarly weighted interval lists " +
            "(by an arbitrary per/base weight), breaking up intervals as needed.") {
    },


    /**
     * A scatter approach that differs from {@link #INTERVAL_SUBDIVISION} in a few ways.
     * <ol>
     * <li>No interval will be subdivided, and consequently, the requested scatter count is an upper bound of scatter count, not a
     * guarantee as to how many {@link IntervalList}s will be produced (e.g., if scatterCount = 10 but there is only one input interval,
     * only 1 interval list will be emitted).</li>
     * <li>When an interval would otherwise be split, it is instead deferred to the next scatter list.</li>
     * <li>The "target width" of each scatter list may be wider than what is computed for INTERVAL_SUBDIVISION.
     * Specifically, if the widest interval in the source interval list is larger than what would otherwise be the target width, that
     * interval's width is used.<br/><br/>The reasoning for this is that this approach produces more consistently-sized interval lists,
     * which is one of the objectives of scattering.</li>
     * </ol>
     */
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION(IntervalListScattererWithoutSubdivision::new, "Scatter the interval list into similarly sized interval lists " +
            "(by base count), but without breaking up intervals."),
    /**
     * A scatter approach that differs from {@link GATKIntervalListScatterMode#BALANCING_WITHOUT_INTERVAL_SUBDIVISION}.
     * <ol>
     * <li>We try to balance the number of unique bases in each interval list by estimating the remaining interval lists sizes.  This is
     * computed from the total number of unique bases and the bases we have consumed.  This means that the interval list with the most
     * number of unique bases is at most the ideal split length larger than the smallest interval list (unique # of bases).</li>
     * </ol>
     */
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW(IntervalListScattererWithoutSubdivisionWithOverflow::new, "Scatter the interval list into similarly sized interval lists " +
            "(by base count), but without breaking up intervals. " +
            "Will overflow current interval list so that the remaining lists will not " +
            "have too many bases to deal with."),

    /**
     * A scatter by interval **count** which attempts to fill each resulting interval list with the same number
     * of intervals, disregarding the base count. This approach can be useful for tools that operate on the interval level
     * rather than the base level, for example CNV calling.
     */
    INTERVAL_COUNT(IntervalListScattererByIntervalCount::new, "Scatter the interval list into similarly sized interval lists " +
            "(by interval count, not by base count). " +
            "Resulting interval lists will contain similar number of intervals.");

    private final Supplier<IntervalListScatterer> scattererSupplier;

    private final String docString;

    GATKIntervalListScatterMode(final Supplier<IntervalListScatterer> supplier, final String docString) {

        scattererSupplier = supplier;
        this.docString = docString;
    }

    @Override
    public String getHelpDoc() {
        return this.docString;
    }

    /**
     * Create the scatterer
     *
     * @return a newly minted Scatterer
     */
    public IntervalListScatterer make() {
        return scattererSupplier.get();
    }
}