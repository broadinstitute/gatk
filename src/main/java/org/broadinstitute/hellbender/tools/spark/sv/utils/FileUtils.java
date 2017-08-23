package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.*;
import java.util.Iterator;

public final class FileUtils {

    public static void writeLinesToSingleFile(final Iterator<String> linesToWrite, final String fileName) {
        try ( final OutputStream writer =
                      new BufferedOutputStream(BucketUtils.createFile(fileName)) ) {
            while (linesToWrite.hasNext()) {
                writer.write(linesToWrite.next().getBytes()); writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write "+fileName, ioe);
        }
    }

    public static void writeSAMFile(final Iterator<SAMRecord> alignments, final SAMFileHeader header, final String outputName,
                                    final boolean preOrdered) {
        try (final SAMFileWriter writer = ReadUtils.createCommonSAMWriter(new File(outputName), null, header,
                preOrdered, true, false) ) {
            alignments.forEachRemaining(writer::addAlignment);
        } catch ( final UncheckedIOException ie) {
            throw new GATKException("Can't write SAM file to the specified location: " + outputName, ie);
        }
    }

    /**
     * Creates a directory, in local FS, HDFS, or Google buckets to write to.
     */
    public static boolean createDirToWriteTo(final String pathString) {
        try {
            Utils.nonNull(pathString);
            if ( java.nio.file.Files.exists(java.nio.file.Paths.get(pathString)) )
                throw new IOException("Directory to be created already exists: " + pathString);

            final boolean isSuccessful;
            if (BucketUtils.isHadoopUrl(pathString)) {
                isSuccessful = org.apache.hadoop.fs.FileSystem.get(new org.apache.hadoop.conf.Configuration()).mkdirs(new org.apache.hadoop.fs.Path(pathString));
            } else {
                final java.nio.file.Path dir = java.nio.file.Files.createDirectory(IOUtils.getPath(pathString));
                isSuccessful = java.nio.file.Files.isDirectory(dir) && java.nio.file.Files.isWritable(dir);
            }
            return isSuccessful;
        } catch (final IOException x) {
            throw new GATKException("Could not create directory: " + pathString, x);
        }
    }
}
