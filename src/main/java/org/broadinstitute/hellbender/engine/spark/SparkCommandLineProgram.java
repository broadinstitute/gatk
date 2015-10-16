package org.broadinstitute.hellbender.engine.spark;

import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.engine.AuthHolder;

import java.io.Serializable;


public abstract class SparkCommandLineProgram extends CommandLineProgram implements Serializable {
    private static final long serialVersionUID = 1l;


    @Argument(doc = "API Key for google cloud authentication",
            shortName = "apiKey", fullName = "apiKey", optional=true)
    protected String apiKey = null;


    @Argument(fullName = "sparkMaster", doc="URL of the Spark Master to submit jobs to when using the Spark pipeline runner.", optional = true)
    protected String sparkMaster = "local[2]";

    @Argument(
            doc = "Name of the program running",
            shortName = "N",
            fullName = "programName",
            optional = true
    )
    public String programName;

    @Override
    protected Object doWork() {
        final JavaSparkContext ctx = SparkContextFactory.getSparkContext(getProgramName(), sparkMaster);
        ctx.getConf().set("spark.driver.maxResultSize", "0")
                .set("spark.driver.userClassPathFirst", "true")
                .set("spark.executor.userClassPathFirst", "true")
                .set("spark.io.compression.codec", "lzf")
                .set("spark.yarn.executor.memoryOverhead", "600");

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

        GCSOptions options = PipelineOptionsFactory.as(GCSOptions.class);
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
}