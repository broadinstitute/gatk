package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.PipelineResult;
import com.google.cloud.dataflow.sdk.options.DataflowPipelineOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.runners.BlockingDataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DirectPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.PipelineRunner;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.Serializable;


public abstract class DataflowCommandLineProgram extends CommandLineProgram implements Serializable {
    private enum PipelineRunnerType implements CommandLineParser.ClpEnum {
        LOCAL(DirectPipelineRunner.class, "run the pipeline locally"),
        BLOCKING(BlockingDataflowPipelineRunner.class, "run the pipeline in the cloud, wait and report status"),
        NONBLOCKING(DataflowPipelineRunner.class, "launch the pipeline in the cloud and don't wait for results");

        public final Class<? extends PipelineRunner<? extends PipelineResult>> runner;
        private final String doc;

        PipelineRunnerType(Class<? extends PipelineRunner<? extends PipelineResult>> runner, String doc){
            this.runner = runner;
            this.doc = doc;
        }

        public String getHelpDoc(){
            return this.doc;
        }

    }
    @Argument(fullName="runner", doc="What pipeline runner to use for dataflow.  Any runner other than LOCAL requires that project and staging be set.")
    private PipelineRunnerType runner = PipelineRunnerType.LOCAL;

    @Argument(fullName="project", doc="dataflow project id", optional=true)
    private String projectID;

    @Argument(fullName = "staging", doc="dataflow staging location, this should be a google bucket of the form gs://", optional = true)
    private String stagingLocation;


    @Override
    protected String[] customCommandLineValidation(){
        if(runner != PipelineRunnerType.LOCAL && (projectID == null || stagingLocation==null)){
            throw new UserException.CommandLineException(String.format("Non local dataflow execution requires project " +
                    "and staging to be set. project:%s id:%s.",projectID, stagingLocation));
        }
        return null;
    }


    @Override
    protected Object doWork() {
        DataflowPipelineOptions options = PipelineOptionsFactory.as(DataflowPipelineOptions.class);
        options.setProject(projectID);
        options.setStagingLocation(stagingLocation);
        options.setRunner(this.runner.runner);
        Pipeline p = Pipeline.create(options);
        DataflowWorkarounds.registerGenomicsCoders(p);
        setupPipeline(p);
        runPipeline(p);
        return null;
    }

    /**
     * Runs a {@link Pipeline} and unwraps RuntimeExceptions caused by UserExceptions back into UserExceptions
     */
    @VisibleForTesting
    protected static void runPipeline(Pipeline p) {
        try{
            p.run();
        } catch( RuntimeException e  ){
            //Data flow catches our UserExceptions and wraps them in RuntimeException
            //Unwrap them again here so that they're handled properly by our help system
            if (e.getCause() instanceof UserException){
                throw (UserException)e.getCause();
            } else {
                throw e;
            }
        }
    }

    /**
     * setup the pipeline for running by adding transforms
     */
    protected abstract void setupPipeline(Pipeline pipeline);
}
