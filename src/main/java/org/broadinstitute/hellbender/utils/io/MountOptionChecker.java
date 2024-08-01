package org.broadinstitute.hellbender.utils.io;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;

public class MountOptionChecker {

    public static boolean pathIsNoExec(final String pathToCheck) throws IOException {
        // We'll check the mount options in /proc/mounts
        try(final BufferedReader mountsReader = new BufferedReader(new FileReader("/proc/mounts"))) {
            String line;
            while((line = mountsReader.readLine()) != null) {
                if(line.contains(pathToCheck)) {
                    String[] splitLine = line.split(" ");
                    if(splitLine.length > 3 && splitLine[1].equals(pathToCheck)) {
                        return splitLine[3].contains("noexec");
                    }
                }
            }
        }
        catch(FileNotFoundException e) {
            // If there's no /proc/mounts file, something weird is probably going on, but it obviously isn't noexec if
            // it doesn't exist, so return false
            return false;
        }
        catch (IOException e) {
            throw e;
        }

        return false;
    }

}
