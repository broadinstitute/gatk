package org.broadinstitute.hellbender.utils;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;

import java.io.Serializable;

/**
 * For command line tools that can implement a Spark and non-Spark version.
 *
 * Subclasses should implement runPipeline(ctx), but ctx can be {@code null}.
 */
public abstract class SparkToggleCommandLineProgram extends SparkCommandLineProgram implements Serializable {
    private static final long serialVersionUID = 1l;
    public static final String DISABLE_SPARK_SHORT_NAME = "ds";
    public static final String DISABLE_SPARK_FULL_NAME = "disableSpark";

    @Argument(
            doc = "Disable spark and run everything in pure Java.",
            shortName = DISABLE_SPARK_SHORT_NAME,
            fullName  = DISABLE_SPARK_FULL_NAME,
            optional  = true
    )
    protected boolean isDisableSpark = false;

    @Override
    protected Object doWork() {

        JavaSparkContext ctx = null;
        if (!isDisableSpark) {
            ctx = SparkContextFactory.getSparkContext(getProgramName(), sparkArgs.getSparkProperties(), sparkArgs.getSparkMaster());
        } else {
            logger.info("Spark disabled.  sparkMaster option (" + sparkArgs.getSparkMaster() + ") ignored.");
        }

        try {
            runPipeline(ctx);
            return null;
        } finally {
            afterPipeline(ctx);
        }
    }

    /**
     * Extend this method to run code after the pipeline returns.
     * This method is called whether or not the runPipeline call succeeded.
     *
     * Please note that the ctx can be null.
     */
    @Override
    protected void afterPipeline(final JavaSparkContext ctx) {
        if (ctx != null) {
            SparkContextFactory.stopSparkContext(ctx);
        }
    }

    /**
     * Returns the program's name.
     * If programName argument is provided, returns that. Otherwise, returns the simple name of the class.
     *
     * Subclasses can override if desired.
     */
    @Override
    protected String getProgramName(){
        return programName == null ? getClass().getSimpleName() : programName;
    }
}
