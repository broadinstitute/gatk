package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.python.StreamingPythonScriptExecutor;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.runtime.AsynchronousStreamWriter;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Example ReadWalker program that uses a Python streaming executor to stream summary data from a BAM
 * input file to a Python process through an asynchronous stream writer. Reads data is accumulated in
 * a List until a batch size threshold is reached, at which point the batch is handed off to the
 * asynchronous stream writer, which writes the batch to the stream on a background thread. The
 * Python process in turn just writes the data to an output file.
 *
 * <ol>
 * <li>Creates a StreamingPythonExecutor.</li>
 * <li>Initializes an AsynchronousWriterService, which is backed by a FIFO file, to allow streaming data to the Python
 * process in batches on a background thread</li>
 * <li>Writes a string of attributes for each read to a List until the batchSize threshold is reached.</li>
 * <li>Starts a background write of the batch, using Python to read each attribute line from the underlying FIFO,
 * and write it back to an output file.</li>
 * </ol>
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
public class ExampleStreamingPythonExecutor extends ReadWalker {
    private final static String NL = System.lineSeparator();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file")
    private File outputFile; // output file produced by Python code

    @Argument(fullName = "batchSize",
            doc = "Size of a batch for writing")
    private int batchSize = 1000;

    // Create the Python executor. This doesn't actually start the Python process, but verifies that
    // the requested Python executable exists and can be located.
    final StreamingPythonScriptExecutor<String> pythonExecutor = new StreamingPythonScriptExecutor<>(true);

    private List<String> batchList = new ArrayList<>(batchSize);
    private int batchCount = 0;

    @Override
    public void onTraversalStart() {

        // Start the Python process, and initialize a stream writer for the executor to use to stream data to Python.
        pythonExecutor.start(Collections.emptyList());
        pythonExecutor.initStreamWriter(AsynchronousStreamWriter.stringSerializer);

        // Execute Python code to open the output file, where it will write the contents of everything it reads
        // from the executor stream. <code sendSynchronousCommand/>
        pythonExecutor.sendSynchronousCommand(String.format("tempFile = open('%s', 'w')" + NL, outputFile.getAbsolutePath()));
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        // Extract data from the read and accumulate, unless we've reached a batch size, in which case we
        // kick off an asynchronous batch write.
        if (batchCount == batchSize) {
            pythonExecutor.waitForPreviousBatchCompletion();
            startAsynchronousBatchWrite();      // start a new batch
        }
        batchList.add(String.format(
                "Read at %s:%d-%d:\n%s\n",
                read.getContig(), read.getStart(), read.getEnd(), read.getBasesString()));
        batchCount++;
    }

    /**
     * On traversal success, write the remaining batch. Post traversal work would be done here.
     * @return Success indicator.
     */
    public Object onTraversalSuccess() {
        pythonExecutor.waitForPreviousBatchCompletion(); // wait for the previous batch to complete, if there is one
        if (batchCount != 0) {
            // If we have any accumulated reads that haven't been dispatched, start one last
            // async batch write, and then wait for it to complete
            startAsynchronousBatchWrite();
            pythonExecutor.waitForPreviousBatchCompletion();
        }

        return true;
    }

    private void startAsynchronousBatchWrite() {
        // Start a new batch write, sending a command to Python to tell it to start consuming the lines that will
        // be written to the stream.
        pythonExecutor.startBatchWrite(
                String.format("for i in range(%s):\n    tempFile.write(tool.readDataFIFO())" + NL + NL, batchCount),
                batchList);
        batchList = new ArrayList<>(batchSize);
        batchCount = 0;
    }

    @Override
    public void closeTool() {
        // Send a synchronous command to Python to close the temp file
        pythonExecutor.sendSynchronousCommand("tempFile.close()" + NL);

        // terminate the executor
        pythonExecutor.terminate();
    }
}
