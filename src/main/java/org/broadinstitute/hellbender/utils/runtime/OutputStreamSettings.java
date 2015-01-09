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
 * Settings that define text to capture from a process stream.
 */
public class OutputStreamSettings {
    private final EnumSet<StreamLocation> streamLocations = EnumSet.noneOf(StreamLocation.class);
    private int bufferSize;
    private File outputFile;
    private boolean appendFile;

    public OutputStreamSettings() {
    }

    /**
     * @param bufferSize The number of bytes to capture, or -1 for unlimited.
     */
    public OutputStreamSettings(int bufferSize) {
        setBufferSize(bufferSize);
    }

    /**
     * @param outputFile The file to write output to.
     */
    public OutputStreamSettings(File outputFile) {
        setOutputFile(outputFile);
    }

    /**
     * @param outputFile The file to write output to.
     * @param append     true if the output file should be appended to.
     */
    public OutputStreamSettings(File outputFile, boolean append) {
        setOutputFile(outputFile, append);
    }

    public OutputStreamSettings(int bufferSize, File outputFile, boolean appendFile) {
        setBufferSize(bufferSize);
        setOutputFile(outputFile, appendFile);
    }

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

    public void clearBufferSize() {
        this.streamLocations.remove(StreamLocation.Buffer);
        this.bufferSize = 0;
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

    public void clearOutputFile() {
        streamLocations.remove(StreamLocation.File);
        this.outputFile = null;
        this.appendFile = false;
    }

    public void printStandard(boolean print) {
        if (print)
            this.streamLocations.add(StreamLocation.Standard);
        else
            this.streamLocations.remove(StreamLocation.Standard);
    }
}
