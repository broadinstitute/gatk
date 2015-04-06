package org.broadinstitute.hellbender.tools.picard.illumina.parser.fakers;

import java.nio.ByteBuffer;

/**
 * Created by jcarey on 3/13/14.
 *
 * Illumina uses an algorithm described in "Theory of RTA" that determines whether or not a cluster passes filter("PF") or not.
 * These values are written as sequential bytes to Filter Files.  The structure of a filter file is as follows:
 * Bytes 0-3  : 0
 * Bytes 4-7  : unsigned int version
 * Bytes 8-11 : unsigned int numClusters
 * Bytes 12...: boolean encoded in a byte representing whether the cluster passing or failing.
 */
public class FilterFileFaker extends FileFaker {

    @Override
    protected void fakeFile(final ByteBuffer buffer) {
        buffer.putInt(0);       // 0
        buffer.putInt(3);       // Version number
        buffer.putInt(1);       // Number of clusters
        buffer.put((byte) 0);   // Byte representing a cluster failing filter(not a PF read)
    }

    @Override
    protected boolean addLeadingZeros() {
        return true;
    }

    @Override
    protected int bufferSize() {
        return Integer.SIZE * 3;
    }
}
