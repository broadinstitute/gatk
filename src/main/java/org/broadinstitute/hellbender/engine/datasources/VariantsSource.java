package org.broadinstitute.hellbender.engine.datasources;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.function.Function;

/**
 * VariantsSource loads variants from a local file (GCS Bucket to come).
 */
public class VariantsSource {

    /**
     * getVariantsList grabs the variants from local files (or perhaps eventually buckets).
     * @param variantSources, list of files  to read from
     */
    public static List<GATKVariant> getVariantsList( List<String> variantSources ) {
        return getVariantsListAs(variantSources, vc -> new VariantContextVariantAdapter(vc));
    }

    /**
     * getVariantsListAs grabs the variants from local files (or perhaps eventually buckets), applies
     * the wrapper function to each object, and returns them as a list of objects of the type returned
     * by the wrapper function
     * @param variantSources list of files  to read from
     * @param wrapFunction function applied to each VariantContext returned
     */
    public static <T> List<T> getVariantsListAs( List<String> variantSources, Function<VariantContext, T> wrapFunction ) {
        final List<T> aggregatedResults = new ArrayList<>();

        for ( final String variantSource : variantSources ) {
            try ( final FeatureDataSource<VariantContext> dataSource =
                          new FeatureDataSource<>(variantSource, null, 0, VariantContext.class) ) {
                aggregatedResults.addAll(wrapQueryResults(dataSource.iterator(), wrapFunction));
            }
        }
        return aggregatedResults;
    }

    private static <T> List<T> wrapQueryResults( final Iterator<VariantContext> queryResults, final Function<VariantContext, T> wrapFunction) {
        final List<T> wrappedResults = new ArrayList<>();
        while ( queryResults.hasNext() ) {
            wrappedResults.add(wrapFunction.apply(queryResults.next()));
        }
        return wrappedResults;
    }

}

