package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.util.stream.StreamSupport;

/**
 * A VariantWalker is a tool that processes a variant at a time from a source of variants, with
 * optional contextual information from a reference and/or sets of reads.
 *
 * VariantWalker authors must implement the apply() method to process each read, and may optionally implement
 * onTraversalStart() and/or onTraversalDone().
 *
 */
public abstract class VariantWalker extends GATKTool {

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF/BCF file containing variants", common = false, optional = false, minElements = 1)
    public File VARIANT_FILE;

    //NOTE: using this directly rather than FeatureInput because this is a special source of variants.
    private FeatureDataSource<VariantContext> variants;

    /**
     * Create the reads and reference data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    @SuppressWarnings("unchecked")
    protected final void onStartup() {
        super.onStartup();

        FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(VARIANT_FILE);
        if (codec.getFeatureType() == VariantContext.class) {
            variants = new FeatureDataSource<>(VARIANT_FILE, (FeatureCodec<VariantContext, ?>)codec);
        } else {
            throw new UserException("File " + VARIANT_FILE + " cannot be decoded as a variant file.");
        }
    }

    /**
     * Implementation of variant-based traversal.
     * Subclasses can override to provide own behavior but default implementation should be suitable for most uses.
     */
    @Override
    public void traverse() {
        VariantFilter filter = makeVariantFilter();
        // Process each variant in the input stream.
        StreamSupport.stream(variants.spliterator(), false)
                .filter(filter)
                .forEach(this::apply);
    }

    /**
     * Returns the variant filter (simple or composite) that will be applied to the variants before calling {@link #apply}.
     * The default implementation filters nothing.
     * Default implementation of {@link #traverse()} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter objects is strongly discouraged.
     *
     * Subclasses can extend to provide own filters (ie override and call super).
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.VariantFilter} composition methods.
     */
    protected VariantFilter makeVariantFilter() {
        return VariantFilterLibrary.ALLOW_ALL_VARIANTS;
    }

    /**
     * Process an individual variant. Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     */
    public abstract void apply( VariantContext variant );

    /**
     * Close the reads and reference data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalDone() instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();

        if ( variants != null )
            variants.close();
    }
}
