package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class VariantsSource {

    private final List<String> variantSources;
    private final Pipeline pipeline;

    public VariantsSource( final List<String> variantSources, final Pipeline pipeline ) {
        for ( final String variantSource : variantSources ) {
            if ( BucketUtils.isCloudStorageUrl(variantSource) ) {
                throw new UnsupportedOperationException("Cloud storage URIs not supported");
            }
        }

        this.variantSources = variantSources;
        this.pipeline = pipeline;
    }

    public PCollection<Variant> getAllVariants() {
        final List<Variant> aggregatedResults = new ArrayList<>();

        for ( final String variantSource : variantSources ) {
            try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(new File(variantSource), getCodecForVariantSource(variantSource), null, 0) ) {
                aggregatedResults.addAll(wrapQueryResults(dataSource.iterator()));
            }
        }

        return pipeline.apply(Create.of(aggregatedResults));
    }


    public PCollection<Variant> getVariantsOverlappingIntervals( final List<SimpleInterval> intervals ) {
        final List<Variant> aggregatedResults = new ArrayList<>();

        for ( final String variantSource : variantSources ) {
            try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(new File(variantSource), getCodecForVariantSource(variantSource), null, 0) ) {
                dataSource.setIntervalsForTraversal(intervals);
                aggregatedResults.addAll(wrapQueryResults(dataSource.iterator()));
            }
        }

        return pipeline.apply(Create.of(aggregatedResults));
    }

    @SuppressWarnings("unchecked")
    private FeatureCodec<VariantContext, ?> getCodecForVariantSource( final String variantSource ) {
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(new File(variantSource));
        if ( !VariantContext.class.isAssignableFrom(codec.getFeatureType()) ) {
            throw new UserException(variantSource + " is not in a format that produces VariantContexts");
        }
        return (FeatureCodec<VariantContext, ?>)codec;
    }

    private List<Variant> wrapQueryResults( final Iterator<VariantContext> queryResults ) {
        final List<Variant> wrappedResults = new ArrayList<>();
        while ( queryResults.hasNext() ) {
            wrappedResults.add(new VariantContextVariantAdapter(queryResults.next()));
        }
        return wrappedResults;
    }
}
