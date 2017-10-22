package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.broadinstitute.hellbender.engine.datasources.ReferenceSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;
import java.util.function.Function;

/**
 * Created by valentin on 10/21/17.
 */
public class ReferenceBackedSVHaplotype extends AbstractSVHaplotype {

    private final ReferenceSource referenceSource;
    private final PipelineOptions pipelineOptions;
    private int length;
   // private final List<AlignmentInterval> alignableIntevals;

    public ReferenceBackedSVHaplotype(final String name, final int length, final List<AlignmentInterval> referenceIntervals,
                                      final ReferenceSource reference, final PipelineOptions pipelineOptions) {
        super(name, referenceIntervals);
        this.length = ParamUtils.isPositive(length, "the length must be positive");
        this.referenceSource = Utils.nonNull(reference);
        this.pipelineOptions = pipelineOptions;
    }


    @Override
    public int getLength() {
        return 0;
    }

    @Override
    public <T> List<List<AlignmentInterval>> align(final Iterable<T> input, final Function<T, byte[]> basesOf) {
        return null;
    }

    @Override
    protected void unsafeCopyBases(int offset, byte[] dest, int destOffset, int length) {

    }
}
