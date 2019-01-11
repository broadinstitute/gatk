package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.CountingVariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.transformers.VariantTransformer;
import org.broadinstitute.hellbender.utils.IndexUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Spliterator;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Base class for variant walkers, which process variants from one or more sources of variants,
 * with optional contextual information from a reference, sets of reads, and/or supplementary sources of
 * Features.
 *
 * Subclasses must implement the {@link #traverse} method to process variants, {@link #initializeDrivingVariants},
 * {@link #getHeaderForVariants},
 * {@link #getSequenceDictionaryForDrivingVariants},
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
     * Returns the pre-filter variant transformer (simple or composite) that will be applied to the variants before filtering.
     * The default implementation uses the {@link VariantTransformer#identity()}.
     * Default implementation of {@link #traverse()} calls this method once before iterating over the variants and reuses
     * the transformer object to avoid object allocation.
     *
     * Subclasses can extend to provide own transformers (i.e. override and call super).
     * Multiple transformers can be composed by using {@link VariantTransformer} composition methods.
     */
    public VariantTransformer makePreVariantFilterTransformer() {
        return VariantTransformer.identity();
    }

    /**
     * Returns the post-filter variant transformer (simple or composite) that will be applied to the variants after filtering.
     * The default implementation uses the {@link VariantTransformer#identity()}.
     * Default implementation of {@link #traverse()} calls this method once before iterating over the variants and reuses
     * the transformer object to avoid object allocation.
     *
     * Subclasses can extend to provide own transformers (i.e. override and call super).
     * Multiple transformers can be composed by using {@link VariantTransformer} composition methods.
     */
    public VariantTransformer makePostVariantFilterTransformer(){
        return VariantTransformer.identity();
    }

    /**
     * Returns a stream over the variants, which are:
     *
     * 1. Transformed with {@link #makePreVariantFilterTransformer()}.
     * 2. Filtered with {@code filter}.
     * 3. Transformed with {@link #makePostVariantFilterTransformer()}.
     */
    protected Stream<VariantContext> getTransformedVariantStream(final CountingVariantFilter filter) {
        final VariantTransformer preTransformer  = makePreVariantFilterTransformer();
        final VariantTransformer postTransformer = makePostVariantFilterTransformer();
        return getTransformedVariantStream(getSpliteratorForDrivingVariants(),
                preTransformer,
                filter,
                postTransformer);
    }

    /**
     * Returns a stream over the variants returned by source, which are:
     *
     * 1. Transformed with preTransformer.
     * 2. Filtered with filter.
     * 3. Transformed with postTransformer.
     */
    protected Stream<VariantContext> getTransformedVariantStream(
            final Spliterator<VariantContext> source,
            final VariantTransformer preTransformer,
            final CountingVariantFilter filter,
            final VariantTransformer postTransformer) {
        return StreamSupport.stream(source, false)
                .map(preTransformer)
                .filter(filter)
                .map(postTransformer);
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
    protected CountingVariantFilter makeVariantFilter() {
        return new CountingVariantFilter(VariantFilterLibrary.ALLOW_ALL_VARIANTS);
    }

}
