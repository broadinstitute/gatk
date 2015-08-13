package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * VariantsDataflowSource loads variants into PCollections using a local file (GCS Bucket to come).
 *
 * This class needs to be mocked, so it cannot be declared final.
 */
public class VariantsDataflowSource {

    private final List<String> variantSources;
    private final Pipeline pipeline;

    /**
     * VariantsDataflowSource sets up source using local files (or eventually GCS buckets).
     * @param variantFiles, list of files (or eventually buckets) to read from
     * @param pipeline, options to get credentials to access GCS buckets.
     */
    public VariantsDataflowSource(final List<String> variantFiles, final Pipeline pipeline) {
        for (final String variantSource : variantFiles) {
            if (BucketUtils.isCloudStorageUrl(variantSource)) {
                // This problem is tracked with issue 632.
                throw new UnsupportedOperationException("Cloud storage URIs not supported");
            }
        }

        this.variantSources = variantFiles;
        this.pipeline = pipeline;
    }

    /**
     * getAllVariants reads variants from file(s) (or eventually GCS bucket(s))
     * @return a PCollection of variants found in the file.
     */
    public PCollection<Variant> getAllVariants() {
        final List<Variant> aggregatedResults = getVariantsList(variantSources);
        if (aggregatedResults.size()==0) {
            // empty list of interval type is something that Dataflow isn't happy with.
            System.err.println("Warning: variant source is empty, you may see a coder failure.");
        }
        return pipeline.apply(Create.of(aggregatedResults)).setName("creatingVariants");
    }

    @SuppressWarnings("unchecked")
    private static FeatureCodec<VariantContext, ?> getCodecForVariantSource( final String variantSource ) {
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(new File(variantSource));
        if ( !VariantContext.class.isAssignableFrom(codec.getFeatureType()) ) {
            throw new UserException(variantSource + " is not in a format that produces VariantContexts");
        }
        return (FeatureCodec<VariantContext, ?>)codec;
    }


    private static List<Variant> wrapQueryResults( final Iterator<VariantContext> queryResults ) {
        final List<Variant> wrappedResults = new ArrayList<>();
        while ( queryResults.hasNext() ) {
            wrappedResults.add(new VariantContextVariantAdapter(queryResults.next()));
        }
        return wrappedResults;
    }


    @VisibleForTesting
    static List<Variant> getVariantsList(List<String> variantSources) {
        final List<Variant> aggregatedResults = new ArrayList<>();

        for ( final String variantSource : variantSources ) {
            try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(new File(variantSource), getCodecForVariantSource(variantSource), null, 0) ) {
                aggregatedResults.addAll(wrapQueryResults(dataSource.iterator()));
            }
        }
        return aggregatedResults;
    }

}

