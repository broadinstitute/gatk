package org.broadinstitute.hellbender.utils;

import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.serializer.KryoSerializer;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;

import java.io.Serializable;

/**
 * For command line tools that can implement a spark and non-spark version.
 *
 * Subclasses should implement runPipeline(ctx), but ctx can be {@code null}.
 */
public abstract class SparkToggleCommandLineProgram extends SparkCommandLineProgram implements Serializable {
    private static final long serialVersionUID = 1l;
    public static final String DISABLE_SPARK_SHORT_NAME = "ds";
    public static final String DISABLE_SPARK_FULL_NAME = "disableSpark";


    @Argument(
            doc = "Skip the spark and run everything in pure Java.",
            shortName = DISABLE_SPARK_SHORT_NAME,
            fullName  = DISABLE_SPARK_FULL_NAME,
            optional  = true
    )
    protected boolean isDisableSpark = false;

    @Override
    protected Object doWork() {

        // TODO: This should be refactored to have the SparkCommandLineProgram share code.  But that would require a change in hellbender as well.
        JavaSparkContext ctx = null;
        if (!isDisableSpark) {
//            ctx = SparkContextFactory.getSparkContext(getProgramName(), sparkMaster);
//            ctx.getConf().set("spark.driver.maxResultSize", "0")
//                    .set("spark.driver.userClassPathFirst", "true")
//                    .set("spark.executor.userClassPathFirst", "true")
//                    .set("spark.io.compression.codec", "lzf")
//                    .set("spark.yarn.executor.memoryOverhead", "600");

            SparkConf sparkConf = new SparkConf().setAppName(getProgramName())
                    .set("spark.serializer", KryoSerializer.class.getCanonicalName())
                    .set("spark.kryo.registrator", "org.broadinstitute.hellbender.engine.spark.GATKRegistrator")
                    .set("spark.kryoserializer.buffer.max", "512m")
                    .set("spark.driver.maxResultSize", "0")
                    .set("spark.driver.userClassPathFirst", "true")
                    .set("spark.executor.userClassPathFirst", "true")
                    .set("spark.io.compression.codec", "lzf")
                    .set("spark.yarn.executor.memoryOverhead", "600")
                    .setMaster(sparkMaster);

            ctx = new JavaSparkContext(sparkConf);
        } else {
            logger.info("Spark disabled.  sparkMaster option (" + this.sparkMaster + ") ignored.");
        }

        try{
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
    protected String getProgramName(){
        return programName == null ? getClass().getSimpleName() : programName;
    }

}
