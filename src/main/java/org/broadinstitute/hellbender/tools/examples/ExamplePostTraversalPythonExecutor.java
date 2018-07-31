package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * Example/toy ReadWalker program that uses a Python script. Processes reads from an input BAM in java,
 * writing read summary data to an output file, then executes a post-traversal Python script process that
 * copies the contents of the java output file in another file.
 *
 * Dependencies (must be included in the GATK docker image):
 * <ul><li>None.</li></ul>
 *
 * See https://github.com/broadinstitute/gatk/wiki/Writing-GATK-Tools-that-use-Python for more information
 * on using Python with GATK.
 */
@CommandLineProgramProperties(
        summary = "Example/toy program that uses a Python script.",
        oneLineSummary = "Example/toy program that uses a Python script.",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public class ExamplePostTraversalPythonExecutor extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file")
    private File javaOutputFile; // output file produced by Java code

    @Argument(fullName = "pythonOutputFile", shortName = "P", doc = "Output file for output of python process")
    private File pythonOutputFile;   // output file produced by Python post-processing step

    // Start the Python executor. This does not actually start the Python process, but fails if python can't be located
    final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);

    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        try {
            outputStream = javaOutputFile != null ? new PrintStream(javaOutputFile) : System.out;
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(javaOutputFile, e);
        }
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        // write summary data to be consumed by the python post-traversal processing step
        outputStream.printf("Read at %s:%d-%d:\n%s\n", read.getContig(), read.getStart(), read.getEnd(), read.getBasesString());
        if ( referenceContext.hasBackingDataSource() )
            outputStream.println("Reference Context:\n" + new String(referenceContext.getBases()));
        outputStream.println();
    }

    @Override
    public Object onTraversalSuccess() {
        // Close the stream, and execute a python script that post-processes the output from the java processing
        // step by just copying it to a separate output file:
        //
        // # Test script that accepts two input arguments that are file paths.
        // # Copies the contents of th first file to the second file:
        // import sys
        // with open(sys.argv[1]) as fin:
        //    lines = fin.readlines()
        //    with open(sys.argv[2], "w") as fout:
        //        fout.writelines(lines)
        //
        outputStream.close();
        outputStream = null;

        final Resource pythonScriptResource = new Resource("copyInputFileToOutputFile.py", ExamplePostTraversalPythonExecutor.class);
        final boolean pythonReturnCode = pythonExecutor.executeScript(
                pythonScriptResource,
                null,
                Arrays.asList(
                        javaOutputFile.getAbsolutePath(),
                        pythonOutputFile.getAbsolutePath()
                )
        );

        return pythonReturnCode;
    }

    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
