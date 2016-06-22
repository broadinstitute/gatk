package org.broadinstitute.hellbender.engine.spark;

import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.scheduler.SchedulerBackend;
import org.apache.spark.scheduler.local.LocalBackend;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.engine.AuthHolder;
import java.io.Serializable;


public abstract class SparkCommandLineProgram extends CommandLineProgram implements Serializable {
    private static final long serialVersionUID = 1L;


    @Argument(doc = "API Key for google cloud authentication",
            shortName = "apiKey", fullName = "apiKey", optional=true)
    protected String apiKey = null;
    JavaSparkContext ctx;

    @Argument(
            doc = "Name of the program running",
            shortName = "N",
            fullName = "programName",
            optional = true
    )
    public String programName;

    @ArgumentCollection
    public SparkCommandLineArgumentCollection sparkArgs = new SparkCommandLineArgumentCollection();



    @Override
    protected void onStartup() {
        super.onStartup();
        ctx = SparkContextFactory.getSparkContext(getProgramName(), sparkArgs.getSparkProperties(), sparkArgs.getSparkMaster());
    }

    @Override
    protected Object doWork() {
        try{
            runPipeline(ctx);
            return null;
        } finally {
            afterPipeline(ctx);
        }
    }

    /**
     * @return a GCSOptions object authenticated with apiKey suitable for accessing files in GCS,
     *         or null if no apiKey is present.
     */
    protected GCSOptions getAuthenticatedGCSOptions() {
        if ( apiKey == null ) {
            return null;
        }

        final GCSOptions options = PipelineOptionsFactory.as(GCSOptions.class);
        options.setApiKey(apiKey);
        return options;
    }

    protected AuthHolder getAuthHolder() {
        return new AuthHolder(getClass().getSimpleName(), apiKey);
    }

    // ---------------------------------------------------
    // Functions meant for overriding

    /**
     * Runs the pipeline.
     */
    protected abstract void runPipeline(final JavaSparkContext ctx);

    /**
     * Extend this method to run code after the pipeline returns.
     * This method is called whether or not the runPipeline call succeeded.
     */
    protected void afterPipeline(final JavaSparkContext ctx) {
        SparkContextFactory.stopSparkContext(ctx);
    }

    /**
     * Returns the program's name.
     * If programName argument is provided, returns that. Otherwise, returns the simple name of the class.
     *
     * Subclasses can override if desired.
     */
    protected String getProgramName(){
        return programName == null ? getClass().getSimpleName() : programName;
    }

    /**
     * Returns the number of cores being used by each worker.
     * In local mode it will return the total number of cores allocated to the process.
     */
    protected int getExecutorCores() {
        int cores;
        SchedulerBackend schedulerBackend = ctx.sc().schedulerBackend();
        if (schedulerBackend instanceof LocalBackend) {
            cores = ((LocalBackend) schedulerBackend).totalCores();
        } else {
            cores = ctx.getConf().getInt("spark.executor.cores",1); // Spark defaults to 1 core per executor if no value is set
        }
        return cores;
    }

    /**
     * Returns the total memory allocated to each core (in MB)
     *
     * NOTE: in local mode this returns the ammount of memory available to the primary executor for caching
     * by default this value will be approximately 0.5 of the available memory for one core. Additionally,
     * there is no guarantee that this result will be correct if the user never specifies either 'spark.executor.memory'
     * or 'spark.executor.cores'
     */
    public double getMemoryPerCore() {
        if (ctx.isLocal()) {
            // In local mode, poll the active executor directly for the memory it has
            return ((long)ctx.sc().getExecutorMemoryStatus().valuesIterator().next()._1() / ((double)getExecutorCores()))/(1000000.0);
        }
        return ctx.sc().executorMemory() / ((double)getExecutorCores());
    }

}
