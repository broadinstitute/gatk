package org.broadinstitute.hellbender.utils.io;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;

import java.io.*;
import java.util.Arrays;

/**
 * Provides utility function(s) for checking the mount options on directories.  Currently limited
 * to checking whether a specific directory is mounted as noexec
 */
public class MountOptionChecker {

    /**
     * Uses the linux/unix cli utility "mount" to check if the specified directory is mounted using the noexec flag
     *
     * @param pathToCheck - String representation of the path to check for a noexec flag
     * @return true if {@code pathToCheck} corresponds to a directory that has been mounted with the noexec option,
     *  otherwise false
     * @throws IOException - BufferedReader.readLine() technically throws this exception, but we're reading from a
     *  BufferedReader wrapped around a StringReader, so if it throws an exception, I'd assume something has gone very
     *  wrong
     * @throws GATKException - if an error occurs running mount
     */
    public static boolean pathIsNoExec(final String pathToCheck) throws IOException, GATKException {
        // Run the mount command without arguments to get mount information
        ProcessController processController = ProcessController.getThreadLocal();
        final String[] command = {"mount"};
        final ProcessSettings settings = new ProcessSettings(command);
        settings.getStdoutSettings().setBufferSize(65536);
        final ProcessOutput mountOutput = processController.exec(settings);

        // Check line by line for the specified path until we find the path and check whether it is noexec
        final BufferedReader stdoutReader = new BufferedReader(new StringReader(mountOutput.getStdout().getBufferString()));
        String line;
        while((line = stdoutReader.readLine()) != null) {
            if(line.contains(pathToCheck)) {
                final String[] splitLine = line.split(" ");
                // Each line of mount output starts with "{filesystem} on {path}" so we'll check the third element for
                // the path
                if(splitLine.length > 3 && splitLine[2].equals(pathToCheck)) {
                    // The flags are in a comma (and space on Mac) separated list enclosed in parentheses at the end of the line
                    final String[] flags = line.substring(line.lastIndexOf('(')+1, line.lastIndexOf(')')).split(", ?");
                    if(Arrays.stream(flags).anyMatch("noexec"::equals)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

}
