/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.runtime;

import java.io.File;
import java.util.Collections;
import java.util.EnumSet;
import java.util.Set;

/**
 * Settings that define text to write to the process stdin.
 */
public class InputStreamSettings {
    private final EnumSet<StreamLocation> streamLocations = EnumSet.noneOf(StreamLocation.class);
    private byte[] inputBuffer;
    private File inputFile;

    public InputStreamSettings() {
    }

    /**
     * @param inputBuffer String to write to stdin.
     */
    public InputStreamSettings(String inputBuffer) {
        setInputBuffer(inputBuffer);
    }

    /**
     * @param inputFile File to write to stdin.
     */
    public InputStreamSettings(File inputFile) {
        setInputFile(inputFile);
    }

    /**
     * @param inputBuffer String to write to stdin.
     * @param inputFile   File to write to stdin.
     */
    public InputStreamSettings(byte[] inputBuffer, File inputFile) {
        setInputBuffer(inputBuffer);
        setInputFile(inputFile);
    }

    public Set<StreamLocation> getStreamLocations() {
        return Collections.unmodifiableSet(streamLocations);
    }

    public byte[] getInputBuffer() {
        return inputBuffer;
    }

    public void setInputBuffer(String inputBuffer) {
        if (inputBuffer == null)
            throw new IllegalArgumentException("inputBuffer cannot be null");
        this.streamLocations.add(StreamLocation.Buffer);
        this.inputBuffer = inputBuffer.getBytes();
    }

    public void setInputBuffer(byte[] inputBuffer) {
        if (inputBuffer == null)
            throw new IllegalArgumentException("inputBuffer cannot be null");
        this.streamLocations.add(StreamLocation.Buffer);
        this.inputBuffer = inputBuffer;
    }

    public void clearInputBuffer() {
        this.streamLocations.remove(StreamLocation.Buffer);
        this.inputBuffer = null;
    }

    public File getInputFile() {
        return inputFile;
    }

    public void setInputFile(File inputFile) {
        if (inputFile == null)
            throw new IllegalArgumentException("inputFile cannot be null");
        this.streamLocations.add(StreamLocation.File);
        this.inputFile = inputFile;
    }

    public void clearInputFile() {
        this.streamLocations.remove(StreamLocation.File);
        this.inputFile = null;
    }

    public void setInputStandard(boolean inputStandard) {
        if (inputStandard)
            this.streamLocations.add(StreamLocation.Standard);
        else
            this.streamLocations.remove(StreamLocation.Standard);
    }
}
