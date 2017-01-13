package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.Arrays;

/**
 * This class represents copy ratio (or copy number) posterior (or prior) expectations
 * as reported by implementations of {@link org.broadinstitute.hellbender.tools.coveragemodel.interfaces.CopyRatioExpectationsCalculator}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CopyRatioExpectations implements Serializable {

    private static final long serialVersionUID = 7573175667965873997L;

    private final double[] logCopyRatioMeans;
    private final double[] logCopyRatioVariances;
    private final int numTargets;

    /**
     * Public constructor.
     *
     * @param logCopyRatioMeans an array of log copy ratio (or copy number) means on a range of targets
     * @param logCopyRatioVariances an array of log copy ratio (or copy number) variances on a range of targets
     */
    public CopyRatioExpectations(@Nonnull final double[] logCopyRatioMeans,
                                 @Nonnull final double[] logCopyRatioVariances) {
        this.logCopyRatioMeans = Utils.nonNull(logCopyRatioMeans, "Log copy ratio means array must be non-null");
        this.logCopyRatioVariances = Utils.nonNull(logCopyRatioVariances, "Log copy ratio variances array must be non-null");
        Utils.validateArg(logCopyRatioMeans.length == logCopyRatioVariances.length, "The length of log copy ratio means and variances arrays must be equal");
        numTargets = logCopyRatioMeans.length;
    }

    /**
     * Return log copy ratio (or copy number) means
     *
     * @return a double array
     */
    public double[] getLogCopyRatioMeans() {
        return logCopyRatioMeans;
    }

    /**
     * Return log copy ratio (or copy number) means on a given target-space block
     *
     * @return a double array
     */
    public double[] getLogCopyRatioMeans(@Nonnull final LinearSpaceBlock block) {
        validateBlockIndexRange(block);
        return Arrays.copyOfRange(logCopyRatioMeans, block.getBegIndex(), block.getEndIndex());
    }

    /**
     * Return log copy ratio (or copy number) variances
     *
     * @return a double array
     */
    public double[] getLogCopyRatioVariances() {
        return logCopyRatioVariances;
    }

    /**
     * Return log copy ratio (or copy number) variances on a given target-space block
     *
     * @return a double array
     */
    public double[] getLogCopyRatioVariances(final LinearSpaceBlock block) {
        validateBlockIndexRange(block);
        return Arrays.copyOfRange(logCopyRatioVariances, block.getBegIndex(), block.getEndIndex());
    }

    private void validateBlockIndexRange(final LinearSpaceBlock block) {
        Utils.nonNull(block, "The target-space block must be non-null");
        ParamUtils.inRange(block.getBegIndex(), 0, numTargets, "Begin index out of range");
        ParamUtils.inRange(block.getEndIndex(), 0, numTargets, "End index out of range");
    }

}
