package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.engine.filters.CountingVariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBOptions;
import org.broadinstitute.hellbender.transformers.VariantTransformer;
import org.broadinstitute.hellbender.utils.IndexUtils;

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
public abstract class VariantWalkerBase extends WalkerBase {

    /**
     * Default value to control the size of the cache for our driving variants input(s)
     * (specifically, the number of additional bases worth of overlapping records to cache for
     * queries on the driving variants).
     */
    public static final int DEFAULT_DRIVING_VARIANTS_LOOKAHEAD_BASES = 100_000;

    //Various options for reading from a GenomicsDB
    protected GenomicsDBOptions genomicsDBOptions;

    @Override
    public boolean requiresFeatures() { return true; }

    @Override
    public String getProgressMeterRecordLabel() { return "variants"; }

    @Override
    protected GenomicsDBOptions getGenomicsDBOptions() {
        if (genomicsDBOptions == null) {
            genomicsDBOptions = new GenomicsDBOptions(referenceArguments.getReferencePath());
        }
        return genomicsDBOptions;
    }

    @Override
    void initializeFeatures() {

        // We override this method to prevent our FeatureManager from being set to null when no side feature inputs
        // are specified, since VariantWalkers always have at least one (driving) variants input. Note that the query
        // lookahead size used here for side inputs is not necessarily the same as for driving variants, which is
        // determined by {@link #getDrivingVariantCacheLookAheadBases}.
        //
        // TODO: I think reducing the lookahead for side inputs from DEFAULT_DRIVING_VARIANTS_LOOKAHEAD_BASES to
        // TODO: FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES will likely hurt performance for tools like VQSR,
        // TODO: but let's test it
        features = new FeatureManager(this, DEFAULT_DRIVING_VARIANTS_LOOKAHEAD_BASES, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                      getGenomicsDBOptions());
        initializeDrivingVariants();
    }

    /**
     * Overriding the superclass method to preferentially
     * choose the sequence dictionary from the driving source of variants.
     */
    @Override
    public SAMSequenceDictionary getBestAvailableSequenceDictionary() {
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
     * When performing a query on the driving variants input(s), the number of additional bases beyond the end
     * of the query for which overlapping variants should be pre-fetched and cached.
     *
     * Defaults to {@link #DEFAULT_DRIVING_VARIANTS_LOOKAHEAD_BASES}
     *
     * Subclasses can customize this value by overriding this method.
     *
     * @return the number of additional bases to prefetch and cache beyond a query on the driving variants
     */
    protected int getDrivingVariantCacheLookAheadBases(){
        return DEFAULT_DRIVING_VARIANTS_LOOKAHEAD_BASES;
    }

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
