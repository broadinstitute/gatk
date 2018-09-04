package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.PicardNonZeroExitException;

/**
 * Adapter shim for use within GATK to run Picard tools. Note that this class does not have
 * it's own CommandLineProgramProperties annotation, and isn't intended to be run from the
 * command line.
 */
public class PicardCommandLineProgramExecutor extends CommandLineProgram {

    // Our wrapped Picard command line program, to which we forward subsequent calls.
    final private picard.cmdline.CommandLineProgram picardCommandLineProgram;

    public PicardCommandLineProgramExecutor(final picard.cmdline.CommandLineProgram picardCommandLineProgram) {
        this.picardCommandLineProgram = picardCommandLineProgram;
    }

    /**
     * Entry point for Picard tools that are called from GATK.
     */
    @Override
    public Object instanceMain(final String[] argv) {
        final int toolReturnCode = picardCommandLineProgram.instanceMain(argv);
        if (toolReturnCode != 0) {
            throw new PicardNonZeroExitException(toolReturnCode);
        }
        return toolReturnCode;
    }

    @Override
    protected Object doWork() {
        // This method should never be called directly. Call instanceMain instead.
        throw new GATKException.ShouldNeverReachHereException(
                String.format("Attempt to call the doWork method on the Picard tool \"%s\" directly.",
                        picardCommandLineProgram.getClass().getName()));
    }
}
