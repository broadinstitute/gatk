package org.broadinstitute.hellbender.engine.spark;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A cache of known sites by file path, with the property that there is only one copy of each collection of known sites per JVM.
 * This class is an alternative for cases that can't use a Spark broadcast due to its 2GB limitation.
 */
class KnownSitesCache {

    private static final Logger log = LogManager.getLogger(KnownSitesCache.class);

    private static final Map<List<String>, IntervalsSkipList<GATKVariant>> PATHS_TO_VARIANTS = new HashMap<>();

    public static synchronized IntervalsSkipList<GATKVariant> getVariants(List<String> paths) {
        if (PATHS_TO_VARIANTS.containsKey(paths)) {
            return PATHS_TO_VARIANTS.get(paths);
        }
        IntervalsSkipList<GATKVariant> variants = retrieveVariants(paths);
        PATHS_TO_VARIANTS.put(paths, variants);
        return variants;
    }

    private static IntervalsSkipList<GATKVariant> retrieveVariants(List<String> paths) {
        return new IntervalsSkipList<>(paths
                .stream()
                .map(KnownSitesCache::loadFromFeatureDataSource)
                .flatMap(Collection::stream)
                .collect(Collectors.toList()));
    }

    private static List<GATKVariant> loadFromFeatureDataSource(String path) {
        int cloudPrefetchBuffer = 40; // only used for GCS
        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(path, null, 0, null, cloudPrefetchBuffer, cloudPrefetchBuffer) ) {
            return wrapQueryResults(dataSource.iterator());
        }
    }

    private static List<GATKVariant> wrapQueryResults(final Iterator<VariantContext> queryResults ) {
        final List<GATKVariant> wrappedResults = new ArrayList<>();
        long count = 0;
        while ( queryResults.hasNext() ) {
            if (count++ % 100000 == 0) {
                log.info("Number of variants read: " + count);
            }
            wrappedResults.add(VariantContextVariantAdapter.sparkVariantAdapter(queryResults.next()));
        }
        return wrappedResults;
    }
}
