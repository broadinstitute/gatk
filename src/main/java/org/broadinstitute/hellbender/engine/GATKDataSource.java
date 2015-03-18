package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Iterator;

/**
 * A GATKDataSource is something that can be iterated over from start to finish
 * and/or queried by genomic interval. It is not necessarily file-based.
 *
 * @param <T> Type of data in the data source
 */
public interface GATKDataSource<T> extends Iterable<T> {
    Iterator<T> query(final SimpleInterval interval);
}
