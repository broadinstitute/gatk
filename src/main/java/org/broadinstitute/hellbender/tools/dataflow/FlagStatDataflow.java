package org.broadinstitute.hellbender.tools.dataflow;


import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.transforms.Combine.AccumulatingCombineFn.Accumulator;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.ReadBAMTransform;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.Contig;
import com.google.cloud.genomics.utils.GenomicsFactory.OfflineAuth;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.FlagStat;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Collectors;

@CommandLineProgramProperties(usage="runs FlagStat on dataflow", usageShort = "FlagStat", programGroup = DataFlowProgramGroup.class)
public class FlagStatDataflow extends DataflowTool {
    private static final Logger LOG = Logger.getLogger(FlagStatDataflow.class.getName());

    @ArgumentCollection
    private IntervalArgumentCollection intervalArgs = new IntervalArgumentCollection();

    @ArgumentCollection
    private RequiredReadInputArgumentCollection readsArguments = new RequiredReadInputArgumentCollection();

   // @ArgumentCollection
   // private RequiredReferenceInputArgumentCollection refArgs = new RequiredReferenceInputArgumentCollection();

    @Argument
    private String outputFile;



    @Override
    protected void setupPipeline(Pipeline pipeline) {
        OfflineAuth auth;
        GCSOptions ops = PipelineOptionsFactory.fromArgs(getCommandLineParser().getArgv()).as(GCSOptions.class);
        GenomicsOptions.Methods.validateOptions(ops);
        try {
            auth = GCSOptions.Methods.createGCSAuth(ops);
        } catch (IOException e) {
            throw new GATKException("a dataflow options issue", e);
        }

        Iterable<Contig> contigs = intervals.stream()
                .map(i -> new Contig(i.getContig(), i.getStart(), i.getEnd()))
                .collect(Collectors.toList());

        List<String> bams = readsArguments.readFiles.stream().map(File::getPath).collect(Collectors.toList());

        PCollection<Read> preads= ReadBAMTransform.getReadsFromBAMFilesSharded(pipeline, auth, contigs, bams);

        //preads.apply(ParDo.of(new readsToCount()));
        preads.apply(Combine.globally(new CombineCounts())).apply(TextIO.Write.to(outputFile));

    }

    private class CombineCounts extends Combine.AccumulatingCombineFn<Read, StatCounter, String> {
        @Override
        public StatCounter createAccumulator() {
            return new StatCounter();
        }

    }


    public class StatCounter implements Accumulator<Read, StatCounter, String> {
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
        public String extractOutput() {
            return stats.toString();
        }
    }

    class readsToCount extends DoFn<Read, FlagStat.FlagStatus> {

        @Override
        public void processElement(ProcessContext processContext) throws Exception {
            FlagStat.FlagStatus stat = new FlagStat.FlagStatus();
            stat.add(processContext.element());

        }

    }





}