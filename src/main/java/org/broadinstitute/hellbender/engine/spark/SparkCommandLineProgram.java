package org.broadinstitute.hellbender.engine.spark;

import org.apache.logging.log4j.Level;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;

import java.io.Serializable;


public abstract class SparkCommandLineProgram extends CommandLineProgram implements Serializable {
    private static final long serialVersionUID = 1L;
    public static final String SPARK_PROGRAM_NAME_LONG_NAME = "program-name";

    @Argument(
            doc = "Name of the program running",
            fullName = SPARK_PROGRAM_NAME_LONG_NAME,
            optional = true
    )
    public String programName;

    @ArgumentCollection
    public SparkCommandLineArgumentCollection sparkArgs = new SparkCommandLineArgumentCollection();


    @Override
    protected Object doWork() {
        final JavaSparkContext ctx = SparkContextFactory.getSparkContext(getProgramName(), sparkArgs.getSparkProperties(), sparkArgs.getSparkMaster());

        // Log the spark arguments for the config here:
        logSparkConfiguration(ctx);

        try{
            runPipeline(ctx);
            return null;
        } finally {
            afterPipeline(ctx);
        }
    }

    // ---------------------------------------------------
    // Private utility methods:
    private void logSparkConfiguration(final JavaSparkContext ctx) {

        // Only log if we would log a debug message:
        if ( logger.getLevel().isMoreSpecificThan(Level.DEBUG) ) {
            logger.debug("Spark Configuration:");
            // I apologize for this use of scala.
            // Apparently Spark is all about the scala and this seems unavoidable...
            for ( final scala.Tuple2<String,String> confArg : ctx.getConf().getAll() ) {
                logger.debug("  " + confArg._1() + " = " + confArg._2());
            }
        }
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
     * If {@link #programName} argument is provided, returns that. Otherwise, returns the simple name of the class.
     *
     * Subclasses can override if desired.
     */
    protected String getProgramName(){
        return programName == null ? getClass().getSimpleName() : programName;
    }
}
