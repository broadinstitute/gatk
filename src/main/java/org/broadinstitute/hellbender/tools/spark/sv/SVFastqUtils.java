package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.apache.hadoop.fs.ContentSummary;
import org.apache.hadoop.fs.FileSystem;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
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
    public static byte[] readToFastqRecord( final GATKRead read, final boolean includeMappingLocation ) {
        final String nameSuffix = read.isPaired() ? (read.isFirstOfPair() ? "/1" : "/2") : "";
        String mapLoc = "";
        if ( includeMappingLocation ) {
            if ( read.isUnmapped() ) mapLoc = " mapping=unmapped";
            else mapLoc = " mapping=" + read.getContig() + ":" + read.getStart();
        }
        final String rec =
                "@" + read.getName() + nameSuffix + mapLoc + "\n" +
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

    /**
     * Load the FASTQ files in the user specified directory and returns an RDD that satisfies the same requirement
     * as described in {@link JavaSparkContext#wholeTextFiles(String, int)}.
     * @param pathToAllInterleavedFASTQFiles path to the directory where all FASTQ files to perform local assembly upon are located
     * @throws GATKException when getting the file count in the specified directory
     */
    public static JavaPairRDD<String, String> loadFASTQFiles(final JavaSparkContext ctx, final String pathToAllInterleavedFASTQFiles){
        try{
            final FileSystem hadoopFileSystem = FileSystem.get(ctx.hadoopConfiguration());
            final ContentSummary cs = hadoopFileSystem.getContentSummary(new org.apache.hadoop.fs.Path(pathToAllInterleavedFASTQFiles));
            final int fileCount = (int) cs.getFileCount();
            return ctx.wholeTextFiles(pathToAllInterleavedFASTQFiles, fileCount);
        }catch (final IOException e){
            throw new GATKException(e.getMessage());
        }
    }
}
