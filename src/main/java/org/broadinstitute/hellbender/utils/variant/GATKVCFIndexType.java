package org.broadinstitute.hellbender.utils.variant;

/**
 * Choose the Tribble indexing strategy
 */
public enum GATKVCFIndexType {
    DYNAMIC_SEEK,       // use DynamicIndexCreator(IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME)
    DYNAMIC_SIZE,       // use DynamicIndexCreator(IndexFactory.IndexBalanceApproach.FOR_SIZE)
    LINEAR,             // use LinearIndexCreator()
    INTERVAL            // use IntervalIndexCreator()
}
