package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval;

import java.io.Closeable;

public interface AnnotatedIntervalWriter extends Closeable {

    /**
     * Write only the header (and any SAMFileHeader or comments)
     *
     * @param annotatedIntervalHeader
     */
    void writeHeader(final AnnotatedIntervalHeader annotatedIntervalHeader);

    /**
     * attempt to close the file
     */
    @Override
    void close();

    /** Write one region to the file.
     *
     * @param annotatedInterval region to write
     */
    void add(final AnnotatedInterval annotatedInterval);
}
