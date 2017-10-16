package org.broadinstitute.hellbender.tools.examples;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.nio.NioBam;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * Example of how to use Spark on Google Cloud Storage directly, without using the GCS Hadoop Connector.
 */
@CommandLineProgramProperties(
    summary = "Example of how to use Spark on Google Cloud Storage directly, without using the GCS Hadoop Connector",
    oneLineSummary = "Example of how to use Spark on Google Cloud Storage directly, without using the GCS Hadoop Connector",
    programGroup = ExampleProgramGroup.class,
    omitFromCommandLine = true
)
public class ExampleNioCountReads extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File OUTPUT_FILE = null;

    @Argument(fullName = "inputPath", shortName = "P", doc = "Input path (eg. gs://foo/bar.bam)", optional = false)
    private String path = null;

    // Typically set to number of executors times number of cores per executor.
    @Argument(fullName = "parts", doc = "number of partitions", optional = false)
    private int parts = 3;

    private void countReads(JavaSparkContext ctx) {
        PrintStream outputStream;

        try {
            outputStream = OUTPUT_FILE != null ? new PrintStream(OUTPUT_FILE) : System.out;
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(OUTPUT_FILE, e);
        }

        NioBam input = new NioBam(path, path + ".bai");
        long readCount = input.getReads(ctx, parts).count();
        outputStream.println("Number of reads: " + readCount);
    }

    /**
     * Runs the pipeline.
     *
     * @param ctx
     */
    @Override
    protected void runPipeline(JavaSparkContext ctx) {
        countReads(ctx);
    }
}
