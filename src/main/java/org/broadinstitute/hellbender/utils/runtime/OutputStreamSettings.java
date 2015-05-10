package org.broadinstitute.hellbender.utils.runtime;

import java.io.File;
import java.util.Collections;
import java.util.EnumSet;
import java.util.Set;

/**
 * Settings that define text to capture from a process stream.
 */
public final class OutputStreamSettings {
    private final EnumSet<StreamLocation> streamLocations = EnumSet.noneOf(StreamLocation.class);
    private int bufferSize;
    private File outputFile;
    private boolean appendFile;

    public Set<StreamLocation> getStreamLocations() {
        return Collections.unmodifiableSet(streamLocations);
    }

    public int getBufferSize() {
        return bufferSize;
    }

    public void setBufferSize(int bufferSize) {
        this.streamLocations.add(StreamLocation.Buffer);
        this.bufferSize = bufferSize;
    }

    public File getOutputFile() {
        return outputFile;
    }

    public boolean isAppendFile() {
        return appendFile;
    }

    /**
     * Overwrites the outputFile with the process output.
     *
     * @param outputFile File to overwrite.
     */
    public void setOutputFile(File outputFile) {
        setOutputFile(outputFile, false);
    }

    public void setOutputFile(File outputFile, boolean append) {
        if (outputFile == null)
            throw new IllegalArgumentException("outputFile cannot be null");
        streamLocations.add(StreamLocation.File);
        this.outputFile = outputFile;
        this.appendFile = append;
    }

    public void printStandard(boolean print) {
        if (print)
            this.streamLocations.add(StreamLocation.Standard);
        else
            this.streamLocations.remove(StreamLocation.Standard);
    }
}
