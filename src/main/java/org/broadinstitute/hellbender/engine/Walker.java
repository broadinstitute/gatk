package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * An abstract class to represent built-in walkers that inherit directly from {@link GATKTool}.
 * Created by jonn on 6/28/18.
 */
public abstract class Walker extends GATKTool {

    @Override
    final protected ReferenceDataSource getReferenceDataSource() {
        throw new GATKException("Should never access ReferenceDataSource in child classes of AssemblyRegionWalker.");
    }

    @Override
    final protected ReadsDataSource getReadsDataSource() {
        throw new GATKException("Should never access ReadsDataSource in child classes of AssemblyRegionWalker.");
    }

    @Override
    final protected FeatureManager getFeatureManager() {
        throw new GATKException("Should never access FeatureManager in child classes of AssemblyRegionWalker.");
    }


}
