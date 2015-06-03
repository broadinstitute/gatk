package org.broadinstitute.hellbender.engine.dataflow;

import com.cloudera.dataflow.spark.SparkPipelineOptions;
import com.cloudera.dataflow.spark.SparkPipelineRunner;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.PipelineResult;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.runners.BlockingDataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DirectPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.PipelineRunner;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;

import java.io.File;
import java.io.Serializable;


public abstract class DataflowCommandLineProgram extends CommandLineProgram implements Serializable {
    private static final long serialVersionUID = 1l;

    protected enum PipelineRunnerType implements CommandLineParser.ClpEnum {
        LOCAL(DirectPipelineRunner.class, "run the pipeline locally"),
        BLOCKING(BlockingDataflowPipelineRunner.class, "run the pipeline in the cloud, wait and report status"),
        NONBLOCKING(DataflowPipelineRunner.class, "launch the pipeline in the cloud and don't wait for results"),
        SPARK(SparkPipelineRunner.class, "run the pipeline on Spark");

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
    @Argument(fullName="runner", doc="What type of pipeline runner to use for dataflow.  " +
        "BLOCKING or NONBLOCKING requires that project and staging be set.")
    protected PipelineRunnerType runnerType = PipelineRunnerType.LOCAL;

    /**
     * Converts the unqualified class name of a runner to the name of the corresponding
     * <code>PipelineRunnerType</code> enum.
     */
    @VisibleForTesting
    public static String getRunnerTypeName(String runnerSimpleName) {
      for (PipelineRunnerType type : PipelineRunnerType.values()) {
        if (type.runner.getSimpleName().equals(runnerSimpleName)) {
          return type.name();
        }
      }
      throw new IllegalArgumentException("No runner found with simple name " + runnerSimpleName);
    }

    @Argument(fullName="project", doc="dataflow project id", optional=true)
    private String projectID;

    @Argument(fullName = "staging", doc="dataflow staging location, this should be a google bucket of the form gs://", optional = true)
    protected String stagingLocation;

    @Argument(doc = "path to the client secrets file for google cloud authentication",
            shortName = "secret", fullName = "client_secret", optional=true, mutex={"apiKey"})
    protected File clientSecret;

    @Argument(doc = "API Key for google cloud authentication",
            shortName = "apiKey", fullName = "apiKey", optional=true, mutex={"client_secret"})
    protected String apiKey = null;

    @Argument(doc = "Number of Dataflow workers to use (or auto if unset).",
            shortName = "numWorkers", fullName = "numWorkers", optional=true)
    protected int numWorkers = 0;

    @Argument(fullName = "sparkMaster", doc="URL of the Spark Master to submit jobs to when using the Spark pipeline runner.", optional = true)
    protected String sparkMaster;

    @Override
    protected String[] customCommandLineValidation(){
        if((runnerType == PipelineRunnerType.BLOCKING || runnerType == PipelineRunnerType.NONBLOCKING)
            && (projectID == null || stagingLocation==null)){
            throw new UserException.CommandLineException(String.format("Non local dataflow execution requires project " +
                    "and staging to be set. project:%s id:%s.",projectID, stagingLocation));
        }
        return null;
    }


    @Override
    protected Object doWork() {
        final Pipeline p = Pipeline.create(buildPipelineOptions());
        DataflowUtils.registerGATKCoders(p);
        setupPipeline(p);
        runPipeline(p);
        afterPipeline(p);
        return null;
    }

    private PipelineOptions buildPipelineOptions() {
        if (sparkMaster == null) {
            // We create GCSOptions instead of DataflowPipelineOptions to keep track of the secrets so we can read
            // data from buckets.
            final GCSOptions options = PipelineOptionsFactory.as(GCSOptions.class);
            options.setProject(projectID);
            options.setStagingLocation(stagingLocation);
            options.setRunner(this.runnerType.runner);
            if (numWorkers!=0) {
                options.setNumWorkers(numWorkers);
            }
            if (apiKey != null) {
                options.setApiKey(apiKey);
            } else if(clientSecret != null) {
                logger.info("Loading " + clientSecret.getName());
                options.setSecretsFile(clientSecret.getAbsolutePath());
            }
            return options;
        } else {
            final SparkPipelineOptions options = PipelineOptionsFactory.as(SparkPipelineOptions.class);
            options.setRunner(this.runnerType.runner);
            options.setSparkMaster(sparkMaster);
            return options;
        }
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


    // ---------------------------------------------------
    // Functions meant for overriding

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
        return runnerType == PipelineRunnerType.BLOCKING || runnerType == PipelineRunnerType.NONBLOCKING;
    }

}
