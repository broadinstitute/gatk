package org.broadinstitute.hellbender.engine.dataflow;

import com.cloudera.dataflow.spark.SparkPipelineOptions;
import com.cloudera.dataflow.spark.SparkPipelineRunner;
import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.PipelineResult;
import com.google.cloud.dataflow.sdk.options.DataflowPipelineOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.runners.BlockingDataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DirectPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.PipelineRunner;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.security.GeneralSecurityException;


public abstract class DataflowCommandLineProgram extends CommandLineProgram implements Serializable {
    private static final long serialVersionUID = 1l;

    // We need authentication options from GCSOptions, and Dataflow options from DataflowPipelineOptions.
    // Neither inherits from the other, so we have to put them together like this.
    //
    // This also includes code to save the OfflineAuth onto the pipelineoptions, so we can get to them later.
    // You see, GCSOptions will store the API Key and the path to client-secrets.json, but the latter
    // doesn't help us if we call createCredentials on the worker (because it won't be able to find
    // the client-secrets file). So here instead we create the OfflineAuth first (which has the
    // contents of the file, if one was specified) and stash that in the options.
    // This guarantees that anyone can call HellbenderDataflowOptions.Methods.getOfflineAuth and
    // succeed.
    public interface HellbenderDataflowOptions extends GCSOptions, DataflowPipelineOptions {

        public static class Methods {
            public static void setOfflineAuth(HellbenderDataflowOptions opts, GenomicsFactory.OfflineAuth auth) throws IOException {
                ByteArrayOutputStream os = new ByteArrayOutputStream();
                try (ObjectOutputStream oos = new ObjectOutputStream(os)) {
                    oos.writeObject(auth);
                    oos.flush();
                }
                opts.setSerializedOfflineAuth(os.toByteArray());
            }
            public static GenomicsFactory.OfflineAuth getOfflineAuth(HellbenderDataflowOptions opts) throws IOException, ClassNotFoundException, GeneralSecurityException {
                byte[] serialized = opts.getSerializedOfflineAuth();
                if (null==serialized && opts.getApiKey()!=null) {
                    // fall back to using the API key only (even if a secrets file was also specified).
                    GenomicsFactory.Builder builder =
                        GenomicsFactory.builder(opts.getAppName()).setNumberOfRetries(opts.getNumberOfRetries());
                    return builder.build().getOfflineAuthFromApiKey(opts.getApiKey());
                }
                try (ObjectInputStream is = new ObjectInputStream(new ByteArrayInputStream(serialized))) {
                    return (GenomicsFactory.OfflineAuth)(is.readObject());
                }
            }

            public static Storage.Objects createStorageClient(HellbenderDataflowOptions opts) throws GeneralSecurityException, IOException, ClassNotFoundException {
                GenomicsFactory.OfflineAuth auth = getOfflineAuth(opts);
                return GCSOptions.Methods.createStorageClient(opts.as(GCSOptions.class), auth);
            }

        }

        // you don't need to call those directly, use the helper methods in Methods instead.
        void setSerializedOfflineAuth(byte[] auth);
        byte[] getSerializedOfflineAuth();

    }

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
            shortName = "apiKey", fullName = "apiKey", optional=true, mutex={"client_secret"}, sensitive = true)
    protected String apiKey = null;

    @Argument(doc = "Number of Dataflow workers to use (or auto if unset).",
            shortName = "numWorkers", fullName = "numWorkers", optional=true)
    protected int numWorkers = 0;

    @Argument(doc = "The Google Compute Engine machine type that Dataflow will use when spinning up worker VMs",
         fullName = "workerMachineType", optional=true)
    protected String workerMachineType = "n1-standard-4";

    @Argument(doc = "Dataflow endpoint to use",
            fullName = "dataflowEndpoint", optional=true)
    protected String dataflowEndpoint = null;

    @Argument(doc = "GCE availability zone for launching Dataflow workers",
            fullName = "zone", optional=true)
    protected String zone = null;


    @Argument(fullName = "sparkMaster", doc="URL of the Spark Master to submit jobs to when using the Spark pipeline runner.", optional = true)
    protected String sparkMaster;

    @Override
    protected String[] customCommandLineValidation() {
        if ((runnerType == PipelineRunnerType.BLOCKING || runnerType == PipelineRunnerType.NONBLOCKING)
            && (projectID == null || stagingLocation==null)) {
            throw new UserException.CommandLineException(String.format("Non local dataflow execution requires project " +
                    "and staging to be set. project:%s staging location:%s.",projectID, stagingLocation));
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

    protected PipelineOptions buildPipelineOptions() {
        if (sparkMaster == null) {
            // We create HellbenderDataflowOptions instead of DataflowPipelineOptions to keep track of the secrets
            // so we can read data from buckets.
            final HellbenderDataflowOptions options = PipelineOptionsFactory.as(HellbenderDataflowOptions.class);
            options.setProject(projectID);
            options.setStagingLocation(stagingLocation);
            options.setRunner(this.runnerType.runner);
            options.setWorkerMachineType(workerMachineType);
            if (null!=dataflowEndpoint) {
                options.setDataflowEndpoint(dataflowEndpoint);
            }
            if (null!=zone) {
                options.setZone(zone);
            }
            // this is new code. If there's a problem, odds are it's our fault and retrying won't help.
            options.setNumberOfRetries(0);
            if (numWorkers!=0) {
                options.setNumWorkers(numWorkers);
            }
            if (apiKey != null) {
                options.setApiKey(apiKey);
            } else if(clientSecret != null) {
                logger.info("Loading " + clientSecret.getName());
                options.setSecretsFile(clientSecret.getAbsolutePath());
            }
            if (apiKey!=null || clientSecret != null) {
                // put a serialized version of the credentials in the pipelineOptions, so we can get to it later.
                try {
                    GenomicsFactory.OfflineAuth offlineAuth = GCSOptions.Methods.createGCSAuth(options);
                    HellbenderDataflowOptions.Methods.setOfflineAuth(options, offlineAuth);
                } catch (Exception x) {
                    throw new GATKException("Error with credentials",x);
                }
            }
            String name = jobName();
            if (null!=name) options.setJobName(name);
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

    /**
     * Override to pick a name for the pipeline.
     * Note that Dataflow requires the name to be unique among running jobs.
     */
    protected String jobName() {
        return null;
    }

    // ---------------------------------------------------
    // Helpers

    /**
     * True if we're configured to run on the cloud.
     */
    public boolean isRemote() {
        return runnerType == PipelineRunnerType.BLOCKING || runnerType == PipelineRunnerType.NONBLOCKING;
    }

}
