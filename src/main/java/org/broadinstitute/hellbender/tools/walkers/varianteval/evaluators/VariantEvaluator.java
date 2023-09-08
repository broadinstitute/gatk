package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

public abstract class VariantEvaluator implements Comparable<VariantEvaluator> {
    private final String simpleName;
    private final VariantEvalEngine engine;

    public VariantEvaluator(VariantEvalEngine engine) {
        this.engine = engine;
        this.simpleName = getClass().getSimpleName();
    }

    //Note: this is used by DISCVR-seq / VariantQC
    public VariantEvaluator(VariantEvalEngine engine, final String simpleName) {
        this.engine = engine;
        this.simpleName = simpleName;
    }

    protected VariantEvalEngine getEngine() {
        return engine;
    }

    // Should return the number of VariantContexts expected as inputs to update.  Can be 1 or 2
    public abstract int getComparisonOrder();

    // called at all sites, regardless of eval context itself; useful for counting processed bases
    // No longer available.  The processed bp is kept in VEW itself for performance reasons
    public void update1(final VariantContext vc, final VariantEvalContext context) {
    }

    public void update2(final VariantContext eval, final VariantContext comp, final VariantEvalContext context) {
    }

    public void finalizeEvaluation() {}

    protected double rate(long n, long d) {
        return n / (1.0 * Math.max(d, 1));
    }

    protected long inverseRate(long n, long d) {
        return n == 0 ? 0 : d / Math.max(n, 1);
    }

    protected double ratio(long num, long denom) {
        return ((double)num) / (Math.max(denom, 1));
    }

    /**
     * Returns true if the variant in vc was a singleton in the original input evaluation
     * set, regardless of variant context subsetting that has occurred.
     * @param eval the VariantContext being assessed for this previous status as a singleton
     * @return true if eval was originally a singleton site
     */
    protected static boolean variantWasSingleton(final VariantContext eval) {
        return eval.getAttributeAsBoolean(VariantEvalEngine.IS_SINGLETON_KEY, false);
    }

    public final String getSimpleName() {
        return simpleName;
    }

    @Override
    public int compareTo(final VariantEvaluator variantEvaluator) {
        return getSimpleName().compareTo(variantEvaluator.getSimpleName());
    }

    /**
     * Evaluation modules that override this function to indicate that they support
     * combining the results of two independent collections of eval data into
     * a single meaningful result.  The purpose of this interface is to
     * allow us to cut up the input data into many independent stratifications, and then
     * at the end of the eval run decide which stratifications to combine.  This is
     * important in the case of AC, where you may have thousands of distinct AC
     * values that chop up the number of variants to too small a number of variants,
     * and you'd like to combine the AC values into ranges containing some percent
     * of the data.
     *
     * For example, suppose you have an eval that
     * counts variants in a variable nVariants.  If you want to be able to combine
     * multiple evaluations of this type, overload the combine function
     * with a function that sets this.nVariants += other.nVariants.
     *
     * Add in the appropriate fields of the VariantEvaluator T
     * (of the same type as this object) to the values of this object.
     *
     * The values in this and other are implicitly independent, so that
     * the values can be added together.
     *
     * @param other a VariantEvaluator of the same type of this object
     */
    public void combine(final VariantEvaluator other) {
        throw new GATKException(getSimpleName() + " doesn't support combining results, sorry");
    }

    /**
     * Must be overloaded to return true for evaluation modules that support the combine operation
     *
     * @return
     */
    public boolean supportsCombine() {
        return false;
    }

    /**
     * Subclasses must overload this to return true if they require an input to include a calling territory, either
     * an interval list or a reference.
     */
    public boolean requiresTerritoryToBeSpecified() {
        return false;
    }
}
