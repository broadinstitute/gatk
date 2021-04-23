package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Base class for pre-packaged walker traversals in the GATK engine.
 *
 * Classes such as {@link ReadWalker} and {@link VariantWalker} in the engine package that implement a standardized traversal should
 * extend this class rather than {@link GATKTool}. Actual concrete tool classes should extend one of the specific walker base
 * types (such as {@link ReadWalker}) or {@link GATKTool} directly rather than this class.
 *
 * Classes that extend {@link WalkerBase} will be unable to directly access engine data sources unless they are in the
 * engine package. This is to allow walker base classes such as {@link ReadWalker} direct datasource access while disallowing
 * it for concrete walker tool implementations, which should get their data via their {@code apply()} method.
 */
public abstract class WalkerBase extends GATKTool {

    /**
     * {@inheritDoc}
     *
     * Walker tools outside of the engine package should not directly access the engine reference datasource.
     * They should get their reference data via {@code apply()} instead. Tools that need direct datasource access
     * (e.g., to implement custom traversal patterns) should extend {@link GATKTool} directly rather than a walker class,
     * or introduce a new walker base class for the new traversal.
     *
     * We are overriding this method to prevent walker tool implementations outside of the engine package from
     * directly accessing the engine datasources, since walker tools should get their data via {@code apply()} instead.
     */
    @Override
    final protected ReferenceDataSource directlyAccessEngineReferenceDataSource() {
        throw new GATKException("Should never directly access the engine ReferenceDataSource in walker tool classes " +
                "outside of the engine package. Walker tools should get their data via apply() instead.");
    }

    /**
     * {@inheritDoc}
     *
     * Walker tools outside of the engine package should not directly access the engine reads datasource.
     * They should get their reads data via {@code apply()} instead. Tools that need direct datasource access
     * (e.g., to implement custom traversal patterns) should extend {@link GATKTool} directly rather than a walker class,
     * or introduce a new walker base class for the new traversal.
     *
     * We are overriding this method to prevent walker tool implementations outside of the engine package from
     * directly accessing the engine datasources, since walker tools should get their data via {@code apply()} instead.
     */
    @Override
    final protected ReadsDataSource directlyAccessEngineReadsDataSource() {
        throw new GATKException("Should never directly access the engine ReadsDataSource in walker tool classes " +
                "outside of the engine package. Walker tools should get their data via apply() instead.");
    }

    /**
     * {@inheritDoc}
     *
     * Walker tools outside of the engine package should not directly access the engine feature manager.
     * They should get their feature data via {@code apply()} instead. Tools that need direct datasource access
     * (e.g., to implement custom traversal patterns) should extend {@link GATKTool} directly rather than a walker class,
     * or introduce a new walker base class for the new traversal.
     *
     * We are overriding this method to prevent walker tool implementations outside of the engine package from
     * directly accessing the engine datasources, since walker tools should get their data via {@code apply()} instead.
     */
    @Override
    final protected FeatureManager directlyAccessEngineFeatureManager() {
        throw new GATKException("Should never directly access the engine FeatureManager in walker tool classes " +
                "outside of the engine package. Walker tools should get their data via apply() instead.");
    }
}
