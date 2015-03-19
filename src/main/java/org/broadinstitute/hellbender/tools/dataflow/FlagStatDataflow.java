package org.broadinstitute.hellbender.tools.dataflow;


import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.CoderRegistry;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.transforms.Combine.AccumulatingCombineFn.Accumulator;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTool;
import org.broadinstitute.hellbender.tools.FlagStat;

import java.util.logging.Logger;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

@CommandLineProgramProperties(usage="runs FlagStat on dataflow", usageShort = "FlagStat", programGroup = DataFlowProgramGroup.class)
public class FlagStatDataflow extends DataflowTool {
    private static final Logger LOG = Logger.getLogger(FlagStatDataflow.class.getName());

    @ArgumentCollection
    private IntervalArgumentCollection intervalArgs = new IntervalArgumentCollection();

    @ArgumentCollection
    private RequiredReadInputArgumentCollection readsArguments = new RequiredReadInputArgumentCollection();

    @ArgumentCollection
    private RequiredReferenceInputArgumentCollection refArgs = new RequiredReferenceInputArgumentCollection();

    @Argument
    private String outputFile;



    @Override
    protected void setupPipeline(Pipeline pipeline) {
        CoderRegistry cr = pipeline.getCoderRegistry();
        cr.registerCoder(SAMRecord.class, SerializableCoder.class);

        ReadsDataSource samReads = new ReadsDataSource(readsArguments.readFiles);
        Stream<Read> reads = StreamSupport.stream(samReads.spliterator(), false)
                .map(ReadConverter::makeRead);


        PCollection< Read > preads = pipeline.apply(Create.of(reads.iterator()));
        //preads.apply(ParDo.of(new readsToCount()));
        preads.apply(Combine.globally(new CombineCounts())).apply(TextIO.Write.to(outputFile));

    }

    private class CombineCounts extends Combine.AccumulatingCombineFn<SAMRecord, StatCounter, String> {
        @Override
        public StatCounter createAccumulator() {
            return new StatCounter();
        }

    }


    public class StatCounter implements Accumulator<SAMRecord, StatCounter, String> {
        private FlagStat.FlagStatus stats = new FlagStat.FlagStatus();

        @Override
        public void addInput(SAMRecord samRecord) {
            stats.add(samRecord);
        }

        @Override
        public void mergeAccumulator(StatCounter statCounter) {
            stats.merge(statCounter.stats);
        }

        @Override
        public String extractOutput() {
            return stats.toString();
        }
    }

    class readsToCount extends DoFn<SAMRecord, FlagStat.FlagStatus> {

        @Override
        public void processElement(ProcessContext processContext) throws Exception {
            FlagStat.FlagStatus stat = new FlagStat.FlagStatus();
            stat.add(processContext.element());

        }

    }





}