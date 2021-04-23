package org.broadinstitute.hellbender;

import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineArgumentValidator;

import java.util.Arrays;

/**
 * Main class to be used as an alternative entry point to {@code org.broadinstitute.hellbender.Main} for performing
 * command line validation only rather than executing the tool. Used only for testing. Failures manifest as exceptions.
 */
public class CommandLineArgumentValidatorMain extends Main {

    // The main entry point to run GATK tools in command line validation mode only.
    public static void main(final String[] argv) {
        new CommandLineArgumentValidatorMain().validateCommandLine(argv);
    }

    /**
     * Call the command line program (specified in the input arguments) in command line validation mode only.
     * @param argv the raw arguments, including the name of the target tool, to run in command line validation mode
     */
    public void validateCommandLine(final String[] argv) {
        final CommandLineProgram program = setupConfigAndExtractProgram(argv, getPackageList(), getClassList(), getCommandLineName());
        final String[] mainArgs = Arrays.copyOfRange(argv, 1, argv.length);
        new CommandLineArgumentValidator(program).instanceMain(mainArgs);
    }
}
