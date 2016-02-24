package org.broadinstitute.hellbender.engine.datasources;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;

import java.io.File;
import java.util.ArrayList;
import java.util.function.Function;
import java.util.Iterator;
import java.util.List;

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
                          new FeatureDataSource<>(new File(variantSource), getCodecForVariantSource(variantSource), null, 0) ) {
                aggregatedResults.addAll(wrapQueryResults(dataSource.iterator(), wrapFunction));
            }
        }
        return aggregatedResults;
    }

    @SuppressWarnings("unchecked")
    private static FeatureCodec<VariantContext, ?> getCodecForVariantSource( final String variantSource ) {
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(new File(variantSource));
        if ( !VariantContext.class.isAssignableFrom(codec.getFeatureType()) ) {
            throw new UserException(variantSource + " is not in a format that produces VariantContexts");
        }
        return (FeatureCodec<VariantContext, ?>)codec;
    }

    private static <T> List<T> wrapQueryResults( final Iterator<VariantContext> queryResults, final Function<VariantContext, T> wrapFunction) {
        final List<T> wrappedResults = new ArrayList<>();
        while ( queryResults.hasNext() ) {
            wrappedResults.add(wrapFunction.apply(queryResults.next()));
        }
        return wrappedResults;
    }

}

