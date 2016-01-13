package org.broadinstitute.hellbender.engine.spark;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.LocalFileSystem;

/**
 * An extension of Hadoop's {@link LocalFileSystem} that doesn't write (or verify) .crc files.
 * This should be used in preference to {@link org.apache.hadoop.fs.RawLocalFileSystem}, since the latter is a not a
 * subclass of {@link LocalFileSystem}, which can cause problems with
 * {@link org.apache.hadoop.fs.FileSystem#getLocal(Configuration)}.
 */
public final class NonChecksumLocalFileSystem extends LocalFileSystem {
    public NonChecksumLocalFileSystem() {
        setWriteChecksum(false);
        setVerifyChecksum(false);
    }
}
