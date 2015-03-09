package org.broadinstitute.hellbender.dataflowtest;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.transforms.Count;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;

import java.io.File;
import java.io.Serializable;

@CommandLineProgramProperties(usage="Count the usages of every word in a text", usageShort="Count Words")
public class WordCount extends CommandLineProgram implements Serializable{
    @Argument
    public File input;

    @Argument
    public File output;

    private void dataflowCount() {

        // Create Pipeline options object.
        PipelineOptions options = PipelineOptionsFactory.create();

        // Create the Pipeline with default options.
        Pipeline p = Pipeline.create(options);

        // Apply a root transform, a text file read, to the pipeline.
        p.apply(TextIO.Read.from(input.getAbsolutePath()))

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
                .apply(TextIO.Write.to(output.getAbsolutePath()));

        p.run();
    }

    @Override
    protected Object doWork() {
        dataflowCount();
        return null;
    }
}
