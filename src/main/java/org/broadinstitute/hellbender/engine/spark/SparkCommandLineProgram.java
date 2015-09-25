package org.broadinstitute.hellbender.engine.spark;

import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;

import java.io.Serializable;


public abstract class SparkCommandLineProgram extends CommandLineProgram implements Serializable {
    private static final long serialVersionUID = 1l;

    @Argument(doc = "API Key for google cloud authentication",
            shortName = "apiKey", fullName = "apiKey", optional=true)
    protected String apiKey = null;


    @Argument(fullName = "sparkMaster", doc="URL of the Spark Master to submit jobs to when using the Spark pipeline runner.", optional = true)
    protected String sparkMaster = "local[2]";

    @Override
    protected Object doWork() {
        final JavaSparkContext ctx = SparkContextFactory.getSparkContext(getProgramName(), sparkMaster);
        ctx.getConf().set("spark.driver.userClassPathFirst", "true")
                     .set("spark.executor.userClassPathFirst", "true")
                     .set("spark.io.compression.codec", "lzf");

        runPipeline(ctx);
        afterPipeline(ctx);

        return null;
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

    // ---------------------------------------------------
    // Functions meant for overriding

    protected abstract void runPipeline(final JavaSparkContext ctx);

    /**
     * Override this to run code after the pipeline returns.
     */
    protected void afterPipeline(final JavaSparkContext ctx) {
        SparkContextFactory.stopSparkContext(ctx);
    }

    protected abstract String getProgramName();
    // ---------------------------------------------------
    // Helpers

}