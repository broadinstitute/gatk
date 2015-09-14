package org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr;

import com.google.cloud.dataflow.sdk.coders.DefaultCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;

import java.util.stream.StreamSupport;

/**
 * Dataflow-ese for "aggregate RecalibrationTables using their 'combine' method".
 * <p>
 * See https://cloud.google.com/dataflow/java-sdk/JavaDoc/com/google/cloud/dataflow/sdk/transforms/Combine.CombineFn
 * for the Dataflow spec this follows. I've copied their javadoc here.
 */
@DefaultCoder(SerializableCoder.class)
public final class RecalibrationTablesMerger extends Combine.CombineFn<RecalibrationTables, RecalibrationTables, RecalibrationTables> {
    private static final long serialVersionUID = 1L;

    public RecalibrationTablesMerger() {
    }

    /**
     * Returns a new, mutable accumulator value, representing the accumulation of zero input values.
     */
    @Override
    public RecalibrationTables createAccumulator() {
        return null;
    }

    /**
     * Adds the given input value to the given accumulator, returning the new accumulator value.
     * For efficiency, the input accumulator may be modified and returned.
     */
    @Override
    public RecalibrationTables addInput(RecalibrationTables accum, RecalibrationTables datum) {
        return combine(accum, datum);
    }

    /**
     * Returns an accumulator representing the accumulation of all the input values accumulated in the merging accumulators.
     * May modify any of the argument accumulators. May return a fresh accumulator, or may return one of the (modified) argument accumulators.
     */
    @Override
    public RecalibrationTables mergeAccumulators(Iterable<RecalibrationTables> accums) {
        return StreamSupport.stream(accums.spliterator(), false)
                .reduce(null, RecalibrationTablesMerger::combine);
    }

    /**
     * Returns the output value that is the result of combining all the input values represented by the given accumulator.
     */
    @Override
    public RecalibrationTables extractOutput(RecalibrationTables tables) {
        return tables;
    }


    private static RecalibrationTables combine(RecalibrationTables accum, RecalibrationTables datum) {
        if (null == accum) {
            return datum;
        } else {
            accum.combine(datum);
            return accum;
        }

    }

}
