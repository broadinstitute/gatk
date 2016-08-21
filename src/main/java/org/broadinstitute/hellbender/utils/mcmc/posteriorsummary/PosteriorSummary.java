package org.broadinstitute.hellbender.utils.mcmc.posteriorsummary;

import java.io.Serializable;

/**
 * Represents central tendency and upper/lower credible-interval bounds of the posterior of a univariate model parameter,
 * along with optional deciles.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PosteriorSummary implements Serializable {

    static final long serialVersionUID = 144L;

    private final double center;
    private final double lower;
    private final double upper;
    private DecileCollection deciles;

    /**
     * Constructs a PosteriorSummary with only given central tendency and upper/lower credible-interval bounds.
     * @param center    central tendency
     * @param lower     lower credible-interval bound
     * @param upper     upper credible-interval bound
     */
    public PosteriorSummary(final double center, final double lower, final double upper) {
        this.center = center;
        this.lower = lower;
        this.upper = upper;
        deciles = null;
    }

    public double getCenter() {   return center; }
    public double getLower() {     return lower;   }
    public double getUpper() {     return upper;   }

    /**
     * Gets the {@link DecileCollection} if it was set previously by {@link PosteriorSummary#setDeciles(DecileCollection)}.
     * @throws IllegalStateException    if deciles were not set previously
     */
    public DecileCollection getDeciles() {
        if (deciles == null) {
            throw new IllegalStateException("Cannot get deciles before they are set.");
        }
        return deciles;
    }

    public void setDeciles(final DecileCollection deciles) {
        this.deciles = deciles;
    }
}
