package org.broadinstitute.hellbender.utils.spark;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.RuntimeIOException;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.read.ReadUtils;


import java.io.*;

/**
 * Miscellaneous Spark-related utilities
 */
public final class SparkUtils {

    private SparkUtils() {}

    /**
     * Converts a headerless Hadoop bam shard (eg., a part0000, part0001, etc. file produced by
     * {@link org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink}) into a readable bam file
     * by adding a header and a BGZF terminator.
     *
     * This method is not intended for use with Hadoop bam shards that already have a header -- these shards are
     * already readable using samtools. Currently {@link ReadsSparkSink} saves the "shards" with a header for the
     * {@link ReadsWriteFormat#SHARDED} case, and without a header for the {@link ReadsWriteFormat#SINGLE} case.
     *
     * @param bamShard The headerless Hadoop bam shard to convert
     * @param header header for the BAM file to be created
     * @param destination path to which to write the new BAM file
     */
    public static void convertHeaderlessHadoopBamShardToBam( final File bamShard, final SAMFileHeader header, final File destination ) {
        try ( FileOutputStream outStream = new FileOutputStream(destination) ) {
            writeBAMHeaderToStream(header, outStream);
            FileUtils.copyFile(bamShard, outStream);
            outStream.write(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
        }
        catch ( IOException e ) {
            throw new UserException("Error writing to " + destination.getAbsolutePath(), e);
        }
    }

    /**
     * Private helper method for {@link #convertHeaderlessHadoopBamShardToBam} that takes a SAMFileHeader and writes it
     * to the provided `OutputStream`, correctly encoded for the BAM format and preceded by the BAM magic bytes.
     *
     * @param samFileHeader SAM header to write
     * @param outputStream stream to write the SAM header to
     */
    private static void writeBAMHeaderToStream( final SAMFileHeader samFileHeader, final OutputStream outputStream ) {
        final BlockCompressedOutputStream blockCompressedOutputStream = new BlockCompressedOutputStream(outputStream, null);
        final BinaryCodec outputBinaryCodec = new BinaryCodec(new DataOutputStream(blockCompressedOutputStream));

        final String headerString;
        final Writer stringWriter = new StringWriter();
        new SAMTextHeaderCodec().encode(stringWriter, samFileHeader, true);
        headerString = stringWriter.toString();

        outputBinaryCodec.writeBytes(ReadUtils.BAM_MAGIC);

        // calculate and write the length of the SAM file header text and the header text
        outputBinaryCodec.writeString(headerString, true, false);

        // write the sequences binarily.  This is redundant with the text header
        outputBinaryCodec.writeInt(samFileHeader.getSequenceDictionary().size());
        for (final SAMSequenceRecord sequenceRecord: samFileHeader.getSequenceDictionary().getSequences()) {
            outputBinaryCodec.writeString(sequenceRecord.getSequenceName(), true, true);
            outputBinaryCodec.writeInt(sequenceRecord.getSequenceLength());
        }

        try {
            blockCompressedOutputStream.flush();
        } catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }
}
