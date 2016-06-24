package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SamPairUtil.PairOrientation;
import org.broadinstitute.hellbender.metrics.MultiLevelMetrics;

import java.io.Serializable;

/**
 * Metrics about the insert size distribution of a paired-end library, created by the
 * CollectInsertSizeMetrics program and usually written to a file with the extension
 * ".insertSizeMetrics".  In addition the insert size distribution is plotted to
 * a file with the extension ".insertSizeMetrics.pdf".
 */
public final class InsertSizeMetrics extends MultiLevelMetrics implements Serializable {

    private static final long serialVersionUID = 1L;

    /**
     *  The MEDIAN insert size of all paired end reads where both ends mapped to the same chromosome.
     */
    public double MEDIAN_INSERT_SIZE;

    /**
     * The median absolute deviation of the distribution.  If the distribution is essentially normal then
     * the standard deviation can be estimated as ~1.4826 * MAD.
     */
    public double MEDIAN_ABSOLUTE_DEVIATION;

    /** The minimum measured insert size.  This is usually 1 and not very useful as it is likely artifactual. */
    public int MIN_INSERT_SIZE;
    /**
     * The maximum measure insert size by alignment. This is usually very high representing either an artifact
     * or possibly the presence of a structural re-arrangement.
     */
    public int MAX_INSERT_SIZE;
    /**
     * The mean insert size of the "core" of the distribution. Artefactual outliers in the distribution often
     * cause calculation of nonsensical mean and stdev values.  To avoid this the distribution is first trimmed
     * to a "core" distribution of +/- N median absolute deviations around the median insert size. By default
     * N=10, but this is configurable.
     */
    public double MEAN_INSERT_SIZE;
    /** Standard deviation of insert sizes over the "core" of the distribution. */
    public double STANDARD_DEVIATION;
    /** The total number of read pairs that were examined in the entire distribution. */
    public long READ_PAIRS;
    /** The pair orientation of the reads in this data category. */
    public PairOrientation PAIR_ORIENTATION;

    /** The "width" of the bins, centered around the median, that encompass 10% of all read pairs. */
    public int WIDTH_OF_10_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 20% of all read pairs. */
    public int WIDTH_OF_20_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 30% of all read pairs. */
    public int WIDTH_OF_30_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 40% of all read pairs. */
    public int WIDTH_OF_40_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 50% of all read pairs. */
    public int WIDTH_OF_50_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 60% of all read pairs. */
    public int WIDTH_OF_60_PERCENT;
    /**
     * The "width" of the bins, centered around the median, that encompass 70% of all read pairs.
     * This metric divided by 2 should approximate the standard deviation when the insert size
     * distribution is a normal distribution.
     */
    public int WIDTH_OF_70_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 80% of all read pairs. */
    public int WIDTH_OF_80_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 90% of all read pairs. */
    public int WIDTH_OF_90_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 100% of all read pairs. */
    public int WIDTH_OF_99_PERCENT;

    /**
     * Return a disambiguating name suffix to be used by the multiple collectors to
     * decorate output names, which are provided by the user in the form of a base name
     * that needs to be disambiguated for each individual collector. This prevents
     * a collector from clobbering the previous collectors output file(s). The value
     * is not used here.
     */
    public static String getUniqueNameSuffix() { return "insertSizeMetrics";}

}
