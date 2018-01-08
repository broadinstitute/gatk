package org.broadinstitute.hellbender.engine.spark;

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
        try{
            runPipeline(ctx);
            return null;
        } finally {
            afterPipeline(ctx);
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
