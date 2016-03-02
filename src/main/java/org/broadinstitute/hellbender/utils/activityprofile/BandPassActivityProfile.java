package org.broadinstitute.hellbender.utils.activityprofile;


import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Set;

/**
 * A band pass filtering version of the activity profile
 *
 * Applies a band pass filter with a Gaussian kernel to the input state probabilities to smooth
 * them out of an interval
 */
public class BandPassActivityProfile extends ActivityProfile {
    public static final int MAX_FILTER_SIZE = 50;
    private final static double MIN_PROB_TO_KEEP_IN_FILTER = 1e-5;
    public static final double DEFAULT_SIGMA = 17.0;

    private final int filterSize;
    private final double sigma;
    private final double[] gaussianKernel;

    /**
     * Create a new BandPassActivityProfile with default sigma and filter sizes
     */
    public BandPassActivityProfile(final Set<SimpleInterval> restrictToIntervals,
                                   final int maxProbPropagationDistance, final double activeProbThreshold,
                                   final SAMFileHeader samHeader) {
        this(restrictToIntervals, maxProbPropagationDistance, activeProbThreshold, MAX_FILTER_SIZE, DEFAULT_SIGMA, samHeader);
    }

    /**
     * sets adaptiveFilterSize to true
     */
    public BandPassActivityProfile(final Set<SimpleInterval> restrictToIntervals,
                                   final int maxProbPropagationDistance, final double activeProbThreshold,
                                   final int maxFilterSize, final double sigma, final SAMFileHeader samHeader) {
        this(restrictToIntervals, maxProbPropagationDistance, activeProbThreshold, maxFilterSize, sigma, true, samHeader);
    }

    /**
     * Create an activity profile that implements a band pass filter on the states
     *
     * @param restrictToIntervals only include states that are within these intervals, if not null
     * @param maxProbPropagationDistance region probability propagation distance beyond it's maximum size
     * @param activeProbThreshold  threshold for the probability of a profile state being active
     * @param maxFilterSize the maximum size of the band pass filter we are allowed to create, regardless of sigma
     * @param sigma the variance of the Gaussian kernel for this band pass filter
     * @param adaptiveFilterSize if true, use the kernel itself to determine the best filter size
     */
    public BandPassActivityProfile(final Set<SimpleInterval> restrictToIntervals, final int maxProbPropagationDistance,
                                   final double activeProbThreshold, final int maxFilterSize, final double sigma, final boolean adaptiveFilterSize,
                                   final SAMFileHeader samHeader) {
        super(maxProbPropagationDistance, activeProbThreshold, restrictToIntervals, samHeader);

        Utils.validateArg(sigma >= 0, "Sigma must be greater than or equal to 0");

        // setup the Gaussian kernel for the band pass filter
        this.sigma = sigma;
        final double[] fullKernel = makeKernel(maxFilterSize, sigma);
        this.filterSize = adaptiveFilterSize ? determineFilterSize(fullKernel, MIN_PROB_TO_KEEP_IN_FILTER) : maxFilterSize;
        this.gaussianKernel = makeKernel(this.filterSize, sigma);
    }

    protected static int determineFilterSize(final double[] kernel, final double minProbToKeepInFilter) {
        final int middle = (kernel.length - 1) / 2;
        int filterEnd = middle;
        while ( filterEnd > 0 ) {
            if ( kernel[filterEnd - 1] < minProbToKeepInFilter ) {
                break;
            }
            filterEnd--;
        }
        return middle - filterEnd;
    }

    protected static double[] makeKernel(final int filterSize, final double sigma) {
        final int bandSize = 2 * filterSize + 1;
        final double[] kernel = new double[bandSize];
        for( int iii = 0; iii < bandSize; iii++ ) {
            kernel[iii] = MathUtils.normalDistribution(filterSize, sigma, iii);
        }
        return MathUtils.normalizeFromRealSpace(kernel);
    }

    /**
     * Our maximize propagation distance is whatever our parent's is, plus our filter size
     *
     * Stops the profile from interpreting sites that aren't yet fully determined due to
     * propagation of the probabilities.
     *
     * @return the distance in bp we might move our probabilities around for some site i
     */
    @Override
    public int getMaxProbPropagationDistance() {
        return super.getMaxProbPropagationDistance() + filterSize;
    }

    /**
     * Get the size (in bp) of the band pass filter
     * @return a positive integer
     */
    public int getBandSize() {
        return 2 * filterSize + 1;
    }

    /**
     * Get the filter size (which is the size of each wing of the band, minus the center point)
     * @return a positive integer
     */
    public int getFilteredSize() {
        return filterSize;
    }

    /**
     * Get the Gaussian kernel sigma value
     * @return a positive double
     */
    public double getSigma() {
        return sigma;
    }

    /**
     * Get the kernel of this band pass filter.  Do not modify returned result
     * @return the kernel used in this band pass filter
     */
    protected double[] getKernel() {
        return gaussianKernel;
    }

    /**
     * Band pass the probabilities in the ActivityProfile, producing a new profile that's band pass filtered
     * @return a new double[] that's the band-pass filtered version of this profile
     */
    @Override
    protected Collection<ActivityProfileState> processState(final ActivityProfileState justAddedState) {
        final Collection<ActivityProfileState> states = new ArrayList<>();

        for ( final ActivityProfileState superState : super.processState(justAddedState) ) {
            if ( superState.isActiveProb() > 0.0 ) {
                for( int i = -filterSize; i <= filterSize; i++ ) {
                    final SimpleInterval loc = getLocForOffset(justAddedState.getLoc(), i);
                    if ( loc != null ) {
                        final double newProb = superState.isActiveProb() * gaussianKernel[i + filterSize];
                        states.add(new ActivityProfileState(loc, newProb));
                    }
                }
            } else {
                states.add(justAddedState);
            }
        }

        return states;
    }
}
