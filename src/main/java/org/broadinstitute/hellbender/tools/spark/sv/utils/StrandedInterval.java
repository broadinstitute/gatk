package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Represents an interval and strand from the reference genome. Users of this class can choose their own interpretation
 * of what strand means. In the context of imprecise variant calling from unassembled breakpoint evidence, strand is set
 * as defined in BreakpointEvidence.isEvidenceUpstreamOfBreakpoint():
 *
 * If true: the evidence suggests a breakpoint at a reference location upstream of the interval's start coordinate
 * If false: the evidence suggests a breakpoint downstream of the interval's end coordinate
 */
@DefaultSerializer(StrandedInterval.Serializer.class)
public class StrandedInterval {
    private static final SVInterval.Serializer intervalSerializer = new SVInterval.Serializer();

    final SVInterval interval;
    final boolean strand;

    public StrandedInterval(final SVInterval interval, final boolean strand) {
        Utils.validate(interval != null, "Can't construct stranded interval with null interval");
        this.interval = interval;
        this.strand = strand;
    }

    public StrandedInterval(final Kryo kryo, final Input input) {
        this.interval = intervalSerializer.read(kryo, input, SVInterval.class);
        this.strand = input.readBoolean();
    }

    private void serialize(final Kryo kryo, final Output output) {
        intervalSerializer.write(kryo, output, interval);
        output.writeBoolean(strand);
    }

    public SVInterval getInterval() {
        return interval;
    }

    public boolean getStrand() {
        return strand;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final StrandedInterval that = (StrandedInterval) o;

        if (strand != that.strand) return false;
        return interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + (strand ? 1 : 0);
        return result;
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<StrandedInterval> {
        @Override
        public void write(final Kryo kryo, final Output output, final StrandedInterval strandedInterval ) {
            strandedInterval.serialize(kryo, output);
        }

        @Override
        public StrandedInterval read(final Kryo kryo, final Input input, final Class<StrandedInterval> klass ) {
            return new StrandedInterval(kryo, input);
        }
    }

    @Override
    public String toString() {
        return "StrandedInterval{" +
                "interval=" + interval +
                ", strand=" + strand +
                '}';
    }
}
