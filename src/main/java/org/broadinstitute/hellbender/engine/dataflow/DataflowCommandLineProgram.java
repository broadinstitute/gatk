package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.PipelineResult;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.runners.BlockingDataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DirectPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.PipelineRunner;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.Serializable;


public abstract class DataflowCommandLineProgram extends CommandLineProgram implements Serializable {
    private static final long serialVersionUID = 1l;

    protected enum PipelineRunnerType implements CommandLineParser.ClpEnum {
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
    @Argument(fullName="runner", doc="What type of pipeline runner to use for dataflow.  Any type other than LOCAL requires that project and staging be set.")
    protected PipelineRunnerType runnerType = PipelineRunnerType.LOCAL;

    @Argument(fullName="project", doc="dataflow project id", optional=true)
    private String projectID;

    @Argument(fullName = "staging", doc="dataflow staging location, this should be a google bucket of the form gs://", optional = true)
    protected String stagingLocation;

    @Argument(doc = "path to the client secrets file for google cloud authentication, necessary if accessing data from buckets.",
            shortName = "secret", fullName = "client_secret", optional=true)
    protected File clientSecret = new File("client-secrets.json");

    @Argument(doc = "API Key for Google Cloud (unnecessary if a client secrets file is provided).",
            shortName = "apiKey", fullName = "apiKey", optional=true)
    protected String apiKey = null;

    @Override
    protected String[] customCommandLineValidation(){
        if(runnerType != PipelineRunnerType.LOCAL && (projectID == null || stagingLocation==null)){
            throw new UserException.CommandLineException(String.format("Non local dataflow execution requires project " +
                    "and staging to be set. project:%s id:%s.",projectID, stagingLocation));
        }
        return null;
    }

    // ---------------------------------------------------
    // Classes meant for overriding

    @Override
    protected Object doWork() {
        // We create GCSOptions instead of DataflowPipelineOptions to keep track of the secrets so we can read
        // data from buckets.
        final GCSOptions options = PipelineOptionsFactory.as(GCSOptions.class);
        options.setProject(projectID);
        options.setStagingLocation(stagingLocation);
        options.setRunner(this.runnerType.runner);
        if (clientSecret!=null) {
            options.setGenomicsSecretsFile(clientSecret.getAbsolutePath()); // TODO: Migrate to secrets file (issue 516).
        }
        if (apiKey!=null) {
            options.setApiKey(apiKey);
        }
        final Pipeline p = Pipeline.create(options);
        DataflowWorkarounds.registerGenomicsCoders(p);
        setupPipeline(p);
        runPipeline(p);
        afterPipeline(p);
        return null;
    }

    /**
     * Runs a {@link Pipeline} and unwraps RuntimeExceptions caused by UserExceptions back into UserExceptions
     */
    @VisibleForTesting
    protected static void runPipeline(final Pipeline p) {
        try{
            p.run();
        } catch( final RuntimeException e  ){
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
     * set up the pipeline for running by adding transforms
     */
    protected abstract void setupPipeline(Pipeline pipeline);

    /**
     * Override this to run code after the pipeline returns.
     */
    protected void afterPipeline(Pipeline pipeline) {}

    // ---------------------------------------------------
    // Helpers

    /**
     * True if we're configured to run on the cloud.
     */
    public boolean isRemote() {
        return runnerType != PipelineRunnerType.LOCAL;
    }

}
