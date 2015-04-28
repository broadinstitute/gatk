package org.broadinstitute.hellbender.tools.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.transforms.Count;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;

@CommandLineProgramProperties(
        usage="Count the usages of every word in a text",
        usageShort="Count Words",
        programGroup = DataFlowProgramGroup.class)
public class WordCount extends DataflowCommandLineProgram {
    @Argument
    public String input;

    @Argument
    public String output;

    @Override
    public void setupPipeline(Pipeline p) {
        // Apply a root transform, a text file read, to the pipeline.
        p.apply(TextIO.Read.from(input))

                // Apply a ParDo transform to the PCollection resulting from the text file read
                .apply(ParDo.of(new DoFn<String, String>() {
                    @Override
                    public void processElement( ProcessContext c ) {
                        String[] words = c.element().split("[^a-zA-Z']+");
                        for ( String word : words ) {
                            if ( !word.isEmpty() ) {
                                c.output(word);
                            }
                        }
                    }
                }))

                        // Apply the Count.PerElement transform to the PCollection of text strings resulting from the ParDo
                .apply(Count.<String>perElement())

                        // Apply a ParDo transform to format the PCollection of word counts from Count() for output
                .apply(ParDo.of(new DoFn<KV<String, Long>, String>() {
                    @Override
                    public void processElement(ProcessContext c) {
                        c.output(c.element().getKey() + ": " + c.element().getValue());
                    }
                }))

                        // Apply a text file write transform to the PCollection of formatted word counts
                .apply(TextIO.Write.to(output));
    }

}
