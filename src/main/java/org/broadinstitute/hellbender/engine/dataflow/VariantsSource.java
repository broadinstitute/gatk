package org.broadinstitute.hellbender.engine.dataflow;

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

    public PCollection<Variant> getVariantCollection( final List<SimpleInterval> intervals ) {
        final List<Variant> variantCollection = new ArrayList<>();

        for ( final String variantSource : variantSources ) {
            variantCollection.addAll(getVariantsForSingleSource(variantSource, intervals));
        }

        return pipeline.apply(Create.of(variantCollection));
    }

    public Iterable<Variant> query(final SimpleInterval interval) {
        return null;
    }
    private List<Variant> getVariantsForSingleSource( final String variantSource, final List<SimpleInterval> intervals ) {
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(new File(variantSource));
        if ( !VariantContext.class.isAssignableFrom(codec.getFeatureType()) ) {
            throw new UserException(variantSource + " is not in a format that produces VariantContexts");
        }

        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(new File(variantSource), (FeatureCodec<VariantContext, ?>) codec) ) {
            dataSource.setIntervalsForTraversal(intervals);

            final List<Variant> variantCollection = new ArrayList<>();
            for ( final VariantContext variant : dataSource ) {
                variantCollection.add(new VariantContextVariantAdapter(variant));
            }

            return variantCollection;
        }
    }
}
