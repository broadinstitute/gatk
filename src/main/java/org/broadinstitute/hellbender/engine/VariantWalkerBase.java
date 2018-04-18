package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.IndexUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Spliterator;
import java.util.function.BiFunction;
import java.util.stream.StreamSupport;

/**
 * Base class for variant walkers, which process one variant at a time from one or more sources of variants,
 * with optional contextual information from a reference, sets of reads, and/or supplementary sources of
 * Features.
 *
 * Subclasses must implement the {@link #apply} method to process each variant, {@link #initializeDrivingVariants},
 * {@link #getHeaderForVariants}, {@link #getSequenceDictionaryForDrivingVariants},
 * {@link #getSpliteratorForDrivingVariants}, and may optionally implement {@link #onTraversalStart},
 * {@link #onTraversalSuccess} and/or {@link #closeTool}.
 */
public abstract class VariantWalkerBase extends GATKTool {

    /**
     * This number controls the size of the cache for our primary and auxiliary FeatureInputs
     * (specifically, the number of additional bases worth of overlapping records to cache when querying feature sources).
     */
    public static final int FEATURE_CACHE_LOOKAHEAD = 100_000;

    @Override
    public boolean requiresFeatures() { return true; }

    @Override
    public String getProgressMeterRecordLabel() { return "variants"; }

    @Override
    void initializeFeatures() {

        //Note: we override this method because we don't want to set feature manager to null if there are no FeatureInputs.
        //This is because we have at least 1 source of features (namely the driving dataset).
        features = new FeatureManager(this, FEATURE_CACHE_LOOKAHEAD, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                      referenceArguments.getReferencePath());
        initializeDrivingVariants();
    }

    /**
     * Overriding the superclass method to preferentially
     * choose the sequence dictionary from the driving source of variants.
     */
    @Override
    public final SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        final SAMSequenceDictionary dictFromDrivingVariants = getSequenceDictionaryForDrivingVariants();
        if (dictFromDrivingVariants != null) {
            //If this dictionary looks like it was synthesized from a feature index, see
            //if there is a better dictionary available from another source (i.e., the reference)
            if (IndexUtils.isSequenceDictionaryFromIndex(dictFromDrivingVariants)) {
                final SAMSequenceDictionary otherDictionary = super.getBestAvailableSequenceDictionary();
                return otherDictionary != null ?
                        otherDictionary :
                        dictFromDrivingVariants;
            } else {
                return dictFromDrivingVariants;
            }
        }
        return super.getBestAvailableSequenceDictionary();
    }

    /**
     * Process the feature inputs that represent the primary driving source(s) of variants for this tool, and
     * perform any necessary header and sequence dictionary validation. Called by the framework during feature
     * initialization.
     */
    protected abstract void initializeDrivingVariants();

    /**
     * Return the VCFHeader to be used for the driving variants for this tool. The value returned will usually
     * have been prepared in {@link #initializeDrivingVariants}
     * @return VCFHeader to be used for the driving variants.
     */
    public abstract VCFHeader getHeaderForVariants();

    /**
     * Return the primary sequence dictionary to be used for the driving variants for this tool. The value returned
     * will usually have been prepared in {@link #initializeDrivingVariants}
     */
    protected abstract SAMSequenceDictionary getSequenceDictionaryForDrivingVariants();

    /**
     * Return a spliterator to be used to iterate over the elements of the driving variants.
     */
    protected abstract Spliterator<VariantContext> getSpliteratorForDrivingVariants();

    /**
     * Implementation of variant-based traversal.
     * Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
     */
    @Override
    public void traverse() {
        final VariantFilter variantfilter = makeVariantFilter();
        final CountingReadFilter readFilter = makeReadFilter();
        // Process each variant in the input stream.
        StreamSupport.stream(getSpliteratorForDrivingVariants(), false)
                .filter(variantfilter)
                .forEach(variant -> {
                    final SimpleInterval variantInterval = new SimpleInterval(variant);
                    apply(variant,
                            new ReadsContext(reads, variantInterval, readFilter),
                            new ReferenceContext(reference, variantInterval),
                            new FeatureContext(features, variantInterval));

                    progressMeter.update(variantInterval);
                });
    }

    /**
     * Returns the variant filter (simple or composite) that will be applied to the variants before calling {@link #apply}.
     * The default implementation filters nothing.
     * Default implementation of {@link #traverse()} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter objects is strongly discouraged.
     *
     * Subclasses can extend to provide own filters (ie override and call super).
     * Multiple filters can be composed by using {@link VariantFilter} composition methods.
     */
    protected VariantFilter makeVariantFilter() {
        return VariantFilterLibrary.ALLOW_ALL_VARIANTS;
    }

    /**
     * Process an individual variant. Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * @param variant Current variant being processed.
     * @param readsContext Reads overlapping the current variant. Will be an empty, but non-null, context object
     *                     if there is no backing source of reads data (in which case all queries on it will return
     *                     an empty array/iterator)
     * @param referenceContext Reference bases spanning the current variant. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current variant's interval
     *                         by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current variant. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply( VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext );

}
