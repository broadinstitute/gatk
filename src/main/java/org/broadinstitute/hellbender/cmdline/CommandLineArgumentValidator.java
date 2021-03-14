package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Adapter shim/alternate GATK entry point for use by GATK tests to run tools in command line argument
 * validation mode. This class does not actually tools, it only validates that the command line arguments
 * are legal for a given invocation.
 *
 * Note that this class does not have it's own CommandLineProgramProperties annotation.
 */
public class CommandLineArgumentValidator extends CommandLineProgram {

    // Our target command line program, to which we delegate arg parsing calls.
    final private CommandLineProgram targetCommandLineProgram;

    public CommandLineArgumentValidator(final CommandLineProgram targetCommandLineProgram) {
        this.targetCommandLineProgram = targetCommandLineProgram;
    }

    /**
     * Entry point to run command line argument validation only.
     */
    @Override
    public Object instanceMain(final String[] argv) {
        if (targetCommandLineProgram instanceof PicardCommandLineProgramExecutor) {
            return ((PicardCommandLineProgramExecutor) targetCommandLineProgram).validateArgs(argv);
        } else {
            // just call parseArgs and then return
            return targetCommandLineProgram.parseArgs(argv);
        }
    }

    @Override
    protected Object doWork() {
        // This method should never be called directly. Call instanceMain instead.
        throw new GATKException.ShouldNeverReachHereException(
                String.format("Attempt to call the doWork method on the validator test tool \"%s\" directly.",
                        targetCommandLineProgram.getClass().getName()));
    }
}
