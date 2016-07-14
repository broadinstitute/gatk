package org.broadinstitute.hellbender.utils.runtime;

import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.Collections;
import java.util.EnumSet;
import java.util.Set;

/**
 * Settings that define text to write to the process stdin.
 */
public final class InputStreamSettings {
    private final EnumSet<StreamLocation> streamLocations = EnumSet.noneOf(StreamLocation.class);
    private byte[] inputBuffer;
    private File inputFile;

    public InputStreamSettings() {
    }

    public Set<StreamLocation> getStreamLocations() {
        return Collections.unmodifiableSet(streamLocations);
    }

    public byte[] getInputBuffer() {
        return inputBuffer;
    }

    public void setInputBuffer(String inputBuffer) {
        Utils.nonNull(inputBuffer,"inputBuffer cannot be null");
        this.streamLocations.add(StreamLocation.Buffer);
        this.inputBuffer = inputBuffer.getBytes();
    }

    public File getInputFile() {
        return inputFile;
    }

    public void setInputFile(File inputFile) {
        Utils.nonNull(inputFile, "inputFile cannot be null");
        this.streamLocations.add(StreamLocation.File);
        this.inputFile = inputFile;
    }
}
