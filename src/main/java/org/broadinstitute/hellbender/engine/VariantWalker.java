package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.Spliterator;

/**
 * A VariantWalker is a tool that processes a variant at a time from a source of variants, with
 * optional contextual information from a reference, sets of reads, and/or supplementary sources
 * of Features.
 *
 * VariantWalker authors must implement the {@link #apply} method to process each variant, and may optionally implement
 * {@link #onTraversalStart}, {@link #onTraversalSuccess} and/or {@link #closeTool}.
 */
public abstract class VariantWalker extends VariantWalkerBase {

    // NOTE: using File rather than FeatureInput<VariantContext> here so that we can keep this driving source
    //       of variants separate from any other potential sources of Features
    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", common = false, optional = false)
    public String drivingVariantFile;

    // NOTE: keeping the driving source of variants separate from other, supplementary FeatureInputs in our FeatureManager in GATKTool
    //we do add the driving source to the Feature manager but we do need to treat it differently and thus this field.
    private FeatureDataSource<VariantContext> drivingVariants;
    private FeatureInput<VariantContext> drivingVariantsFeatureInput;

    @Override
    protected SAMSequenceDictionary getSequenceDictionaryForDrivingVariants() { return drivingVariants.getSequenceDictionary(); }

    @Override
    protected Spliterator<VariantContext> getSpliteratorForDrivingVariants() { return drivingVariants.spliterator(); }

    /**
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalStart} instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();
        if ( hasUserSuppliedIntervals() ) {
            drivingVariants.setIntervalsForTraversal(userIntervals);
        }
    }

    @Override
    protected void initializeDrivingVariants() {
        drivingVariantsFeatureInput = new FeatureInput<>(drivingVariantFile, "drivingVariantFile");

        //This is the data source for the driving source of variants, which uses a cache lookahead of FEATURE_CACHE_LOOKAHEAD
        drivingVariants = new FeatureDataSource<>(drivingVariantsFeatureInput, FEATURE_CACHE_LOOKAHEAD, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                                  referenceArguments.getReferencePath());

        //Add the driving datasource to the feature manager too so that it can be queried. Setting lookahead to 0 to avoid caching.
        //Note: we are disabling lookahead here because of windowed queries that need to "look behind" as well.
        features.addToFeatureSources(0, drivingVariantsFeatureInput, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                     referenceArguments.getReferencePath());

        //Note: the intervals for the driving variants are set in onStartup
    }

    /**
     * Returns the feature input for the driving variants file.
     */
    protected final FeatureInput<VariantContext> getDrivingVariantsFeatureInput() {
        return drivingVariantsFeatureInput;
    }

    /**
     * Gets the header associated with our driving source of variants as a VCFHeader.
     *
     * @return VCFHeader for our driving source of variants
     */
    public final VCFHeader getHeaderForVariants() {
        final Object header = drivingVariants.getHeader();

        if ( ! (header instanceof VCFHeader) ) {
            throw new GATKException("Header for " + drivingVariantFile + " is not in VCF header format");
        }

        return (VCFHeader)header;
    }

    /**
     * Close all data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalSuccess} and/or
     * {@link #closeTool} instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();

        if ( drivingVariants != null )
            drivingVariants.close();
    }
}
