package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.PipelineResult;
import com.google.cloud.dataflow.sdk.options.DataflowPipelineOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.runners.BlockingDataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DirectPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.PipelineRunner;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.Serializable;


public abstract class DataflowTool extends CommandLineProgram implements Serializable {
    private enum PipelineRunnerType {
        LOCAL(DirectPipelineRunner.class),
        BLOCKING(BlockingDataflowPipelineRunner.class),
        NONBLOCKING(DataflowPipelineRunner.class);

        public  Class<? extends PipelineRunner<? extends PipelineResult>> runner;

        private PipelineRunnerType(Class<? extends PipelineRunner<? extends PipelineResult>> runner){
            this.runner = runner;
        }

    }
    @Argument(fullName="runner", doc="What pipeline runner to use for dataflow.  Any runner other than LOCAL requires that project and staging be set.")
    PipelineRunnerType runner = PipelineRunnerType.LOCAL;

    @Argument(fullName="project", doc="dataflow project id", optional=true)
    String projectID;

    @Argument(fullName = "staging", doc="dataflow staging location, this should be a google bucket of the form gs://", optional = true)
    String stagingLocation;


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
        setupPipeline(p);
        p.run();
        return null;
    }

    /**
     * setup the pipeline for running by adding transforms
     */
    protected abstract void setupPipeline(Pipeline pipeline);
}
