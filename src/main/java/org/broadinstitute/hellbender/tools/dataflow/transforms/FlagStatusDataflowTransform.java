package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.tools.FlagStat;

import java.io.Serializable;

/**
 * Computes Flag stats a {@link PCollection<Read>}
 */
public final class FlagStatusDataflowTransform extends PTransformSAM<FlagStat.FlagStatus> {
    @Override
    public PCollection<FlagStat.FlagStatus> apply(PCollection<Read> input) {
        return input.apply(Combine.globally(new CombineCounts()));
    }

    private static class CombineCounts extends Combine.AccumulatingCombineFn<Read, StatCounter, FlagStat.FlagStatus> {
        @Override
        public StatCounter createAccumulator() {
            return new StatCounter();
        }

    }


    private static class StatCounter implements Combine.AccumulatingCombineFn.Accumulator<Read, StatCounter, FlagStat.FlagStatus>, Serializable {
        private FlagStat.FlagStatus stats = new FlagStat.FlagStatus();

        @Override
        public void addInput(Read read) {
            stats.add(read);
        }

        @Override
        public void mergeAccumulator(StatCounter statCounter) {
            stats.merge(statCounter.stats);
        }

        @Override
        public FlagStat.FlagStatus extractOutput() {
            return stats;
        }
    }

}
