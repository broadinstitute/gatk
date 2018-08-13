package org.broadinstitute.hellbender.engine.datasources;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * A cache of reference bases by file path, with the property that there is only one copy of each collection per JVM.
 * This class is an alternative for cases that can't use a Spark broadcast due to its 2GB limitation.
 */
public class ReferenceSourceCache {

    private static final Logger log = LogManager.getLogger(ReferenceSourceCache.class);

    private static final Map<String, ReferenceSource> PATH_TO_REFERENCE_SOURCES = new HashMap<>();

    public static synchronized ReferenceSource getReferenceSource(String path) throws IOException {
        if (PATH_TO_REFERENCE_SOURCES.containsKey(path)) {
            return PATH_TO_REFERENCE_SOURCES.get(path);
        }
        ReferenceSource referenceSource = retrieveReferenceSource(path);
        PATH_TO_REFERENCE_SOURCES.put(path, referenceSource);
        return referenceSource;
    }

    private static ReferenceSource retrieveReferenceSource(String path) throws IOException {
        return new ReferenceCachingSource(path);
    }
}