package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.python.StreamingPythonScriptExecutor;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.runtime.AsynchronousStreamWriterService;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 * Example ReadWalker program that uses a Python streaming executor to stream summary data from a BAM
 * input file to a Python process through a FIFO. The read data is accumulated in a List until a batch
 * size threshold is reached, at which point the batch is handed off to an asynchronous write service,
 * which writes the batch to the FIFO stream on a background thread. The Python process in turn just
 * writes the data to an output file.
 *
 * <ol>
 * <li>Opens a FIFO for writing.</li>
 * <li>Creates an AsynchronousWriterService to allow writing to the FIFO in batches on a background thread</li>
 * <li>Writes a string of attributes for each read to the List until the batchSize threshold is reached.</li>
 * <li>Uses Python to read each attribute line from the FIFO, and write it to the output file.</li>
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
    final StreamingPythonScriptExecutor pythonExecutor = new StreamingPythonScriptExecutor(true);

    private FileOutputStream fifoWriter;
    private AsynchronousStreamWriterService<String> asyncWriter = null;
    private List<String> batchList = new ArrayList<>(batchSize);
    private int batchCount = 0;

    @Override
    public void onTraversalStart() {

        // Start the Python process, and get a FIFO from the executor to use to send data to Python. The lifetime
        // of the FIFO is managed by the executor; the FIFO will be destroyed when the executor is terminated.
        pythonExecutor.start(Collections.emptyList());
        final File fifoFile = pythonExecutor.getFIFOForWrite();

        // Open the FIFO for writing. Opening a FIFO for read or write will block until there is reader/writer
        // on the other end, so before we open it, send an ASYNCHRONOUS command, that doesn't wait for a
        // response, to the Python process to open the FIFO for reading. The Python process will then block until
        // we open the FIFO.
        pythonExecutor.sendAsynchronousCommand(String.format("fifoFile = open('%s', 'r')" + NL, fifoFile.getAbsolutePath()));
        try {
            fifoWriter = new FileOutputStream(fifoFile);
            asyncWriter = pythonExecutor.getAsynchronousStreamWriterService(fifoWriter, AsynchronousStreamWriterService.stringSerializer);
        } catch ( IOException e ) {
            throw new GATKException("Failure opening FIFO for writing", e);
        }
        // synchronize on the output prompt before executing the next statement
        pythonExecutor.getAccumulatedOutput();

        // Also, ask Python to open our output file, where it will write the contents of everything it reads
        // from the FIFO. <code sendSynchronousCommand/>
        pythonExecutor.sendSynchronousCommand(String.format("tempFile = open('%s', 'w')" + NL, outputFile.getAbsolutePath()));
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        // Extract data from the read and accumulate, unless we've reached a batch size, in which case we
        // kick off an asynchronous batch.
        if (batchCount == batchSize) {
            startAsynchronousBatchWrite();
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
        if (batchCount != 0) {
            // we have accumulated reads that haven't been dispatched;, so dispatch the last batch
            startAsynchronousBatchWrite();
        }
        // make sure the final batch has completed
        asyncWriter.waitForPreviousBatchCompletion(1000, TimeUnit.MILLISECONDS);

        // synchronize on the output prompt before executing the next statement
        pythonExecutor.getAccumulatedOutput();

        // Send synchronous commands to Python to close the temp file and the FIFO file
        // Terminate the async writer and Python executor in closeTool, since this always gets called.
        pythonExecutor.sendSynchronousCommand("tempFile.close()" + NL);
        pythonExecutor.sendSynchronousCommand("fifoFile.close()" + NL);

        return true;
    }

    private void startAsynchronousBatchWrite() {
        // Before we hand off a new batch, wait for the previous batch to complete.
        asyncWriter.waitForPreviousBatchCompletion(1000, TimeUnit.MILLISECONDS);

        // Send an ASYNCHRONOUS command to Python to tell it to start consuming the lines about to be written
        // to the FIFO. Sending a *SYNCHRONOUS* command here would immediately block the background thread
        // since this statement will be executed BEFORE any data from the batch is written to the stream.
        pythonExecutor.sendAsynchronousCommand(String.format(
                "for i in range(%s):\n    tempFile.write(fifoFile.readline())" + NL + NL, batchCount));
        asyncWriter.startAsynchronousBatchWrite(batchList);
        batchList = new ArrayList<>(batchSize);
        batchCount = 0;
    }

    @Override
    public void closeTool() {
        if (asyncWriter != null) {
            asyncWriter.terminate();
        }
        try {
            fifoWriter.close();
        } catch (IOException e) {
            throw new GATKException("IOException closing fifo writer", e);
        }
        pythonExecutor.terminate();
    }
}
