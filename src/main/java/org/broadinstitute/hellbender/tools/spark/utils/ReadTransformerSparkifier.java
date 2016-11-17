package org.broadinstitute.hellbender.tools.spark.utils;

import org.apache.spark.api.java.function.Function;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class ReadTransformerSparkifier implements Function<GATKRead, GATKRead> {

    private static final long serialVersionUID = 1L;

    private final ReadTransformer transformer;

    public ReadTransformerSparkifier( final ReadTransformer transformer ) {this.transformer = transformer;}

    @Override public GATKRead call( final GATKRead read ) { return transformer.apply(read); }

}
