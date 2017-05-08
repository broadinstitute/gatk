package org.broadinstitute.hellbender.tools.spark.utils;

import java.util.Collection;

/**
 * Interface for super-sets containing QueryableLongSets
 */
public interface LargeQueryableLongSet extends QueryableLongSet {

    Collection<? extends QueryableLongSet> getSets();

}
