package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

/**
 * Memory-economical utilities for producing a FASTQ file.
 */
public class SVFastqUtils {

    /** Convert a read into a FASTQ record, represented as a byte[]. */
    public static byte[] readToFastqRecord( final GATKRead read ) {
        final String nameSuffix = read.isPaired() ? (read.isFirstOfPair() ? "/1" : "/2") : "";
        final String rec =
                "@" + read.getName() + nameSuffix + "\n" +
                read.getBasesString() + "\n" +
                "+\n" +
                ReadUtils.getBaseQualityString(read)+"\n";
        return rec.getBytes();
    }

    /**
     * Sort a list of FASTQ records.  (Probably ought to be an ArrayList for memory-efficiency.)
     * This puts them into proper order for an interleaved FASTQ file (since the FASTQ record begins with the
     * template name).
     */
    public static void sortFastqRecords( final List<byte[]> fastqRecordList ) {
        fastqRecordList.sort(SVFastqUtils::compareByteArrays);
    }

    /** Write a list of FASTQ records into a file. */
    public static void writeFastqFile(
            final String fileName,
            final PipelineOptions pipelineOptions,
            List<byte[]> fastqRecordList ) {
        try ( final OutputStream writer =
                      new BufferedOutputStream(BucketUtils.createFile(fileName, pipelineOptions)) ) {
            for ( final byte[] fastqRecord : fastqRecordList ) {
                writer.write(fastqRecord);
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write "+fileName, ioe);
        }
    }

    /** Sure seems like this should be available in the standard Java library somewhere, but I can't find it. */
    private static int compareByteArrays( final byte[] arr1, final byte[] arr2 ) {
        final int len = Math.min(arr1.length, arr2.length);
        for ( int idx = 0; idx != len; ++idx ) {
            final int result = Integer.compare(arr1[idx] & 0xff, arr2[idx] & 0xff);
            if ( result != 0 ) return result;
        }
        return Integer.compare(arr1.length, arr2.length);
    }
}
