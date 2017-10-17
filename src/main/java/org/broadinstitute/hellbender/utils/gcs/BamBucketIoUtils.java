package org.broadinstitute.hellbender.utils.gcs;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Utility methods to deal with read-alignment files stored in buckets
 * or local file system.
 *
 * @author &lt;valentin@broadinstitute.org&gt; on 7/24/17.
 */
public class BamBucketIoUtils {

    /**
     * Creates a writer to a new read-aligment file given its destination resource name and
     * header.
     * <p>
     *     The format of the output will depend on the file name. Currently only BAM (*.bam) and SAM (*.sam)
     *     formats are supported.
     * </p>
     * @param dest the destination file name; it might make reference to a bucket location, hadoop file system resource or a
     *             regular file.
     * @param header the output header. It will indicate the expected order of the read-alignment records in
     *                the output.
     * @param preOrdered whether the calling code will add each record in the correct order or whether the
     *                   writer must sort them itself.
     * @return never {@code null}
     */
    public static SAMFileWriter makeWriter(final String dest, final SAMFileHeader header, final boolean preOrdered) {
        Utils.nonNull(dest, "the destination string cannot be null");
        Utils.nonNull(header, "the input header cannot be null");
        // Using the default factory configuration.
        // perhaps in the future we shall allow to configure this:
        final SAMFileWriterFactory factory = new SAMFileWriterFactory();
        final int lastDotIndex = dest.lastIndexOf('.');
        if (lastDotIndex >= 0) {
            final String extension = dest.substring(lastDotIndex).toLowerCase();
            if (extension.equals(BamFileIoUtils.BAM_FILE_EXTENSION)) {
                return factory.makeBAMWriter(header, preOrdered, BucketUtils.createFile(dest));
            } else if (extension.equals(".sam")) {
                return factory.makeSAMWriter(header, preOrdered, BucketUtils.createFile(dest));
            } else {
                throw new GATKException("unsupported read alignment file name extension (." + extension + ") in requested name: " + dest);
            }
        } else {
            throw new GATKException("cannot determine the alignment file format from its name: " + dest);
        }
    }
}
