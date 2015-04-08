package org.broadinstitute.hellbender.tools.dataflow;


import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.Sum;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.DataFlowSAMFn;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;


@CommandLineProgramProperties(usage = "Count bases in dataflow with filters and transformers", usageShort = "count bases", programGroup = DataFlowProgramGroup.class)
public class CountBasesFiltered extends DataflowReadsPipeline {

    public CountBasesFiltered() {
        super(new CountBases(), ImmutableList.of(r -> r.getReferenceName().equals("chr7")), ImmutableList.of());
    }

    public static class CountBases extends PTransformSAM {
        @Override
        public PCollection<String> apply(PCollection<Read> reads) {

            return reads.apply(ParDo.of(new DataFlowSAMFn<Long>(getHeaderString()) {
                @Override
                protected void apply(SAMRecord read) {
                    Long bases = (long) read.getReadBases().length;
                    output(bases);
                }
            }))
                .apply(Sum.longsGlobally())
                .apply(ParDo.of(new DoFn<Long, String>() {
                    @Override
                    public void processElement(DoFn<Long, String>.ProcessContext c) throws Exception {
                        c.output(String.valueOf(c.element()));
                }
            }));
        }

    }
}


