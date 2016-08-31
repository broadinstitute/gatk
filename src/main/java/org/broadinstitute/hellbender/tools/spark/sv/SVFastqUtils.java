package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.hadoop.fs.ContentSummary;
import org.apache.hadoop.fs.FileSystem;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Memory-economical utilities for producing a FASTQ file.
 */
public class SVFastqUtils {

    /** Convert a read into a FASTQ record, represented as a byte[]. */
    @VisibleForTesting
    public static byte[] readToFastqRecord( final GATKRead read, final boolean includeMappingLocation ) {
        final String nameSuffix = read.isPaired() ? (read.isFirstOfPair() ? "/1" : "/2") : "";
        String mapLoc = "";
        if ( includeMappingLocation ) {
            if ( read.isUnmapped() ) mapLoc = " mapping=unmapped";
            else mapLoc = " mapping=" + read.getContig() + ":" + read.getStart() + ";" + read.getCigar().toString();
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
    @VisibleForTesting
    public static void sortFastqRecords( final List<byte[]> fastqRecordList ) {
        fastqRecordList.sort(SVFastqUtils::compareByteArrays);
    }

    /** Write a list of FASTQ records into a file. */
    @VisibleForTesting
    public static void writeFastqFile(final String fileName,
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
    @VisibleForTesting
    static int compareByteArrays( final byte[] arr1, final byte[] arr2 ) {
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
     * @throws IllegalArgumentException if the input directory is empty
     */
    @VisibleForTesting
    public static JavaPairRDD<String, String> loadFASTQFiles(final JavaSparkContext ctx, final String pathToAllInterleavedFASTQFiles){
        try{
            final int fileCount = (int) FileSystem.get(ctx.hadoopConfiguration()).getContentSummary(new org.apache.hadoop.fs.Path(pathToAllInterleavedFASTQFiles)).getFileCount();
            if(fileCount==0){
                throw new IllegalArgumentException("Empty directory: " + pathToAllInterleavedFASTQFiles);
            }
            return ctx.wholeTextFiles(pathToAllInterleavedFASTQFiles, fileCount);
        }catch (final IOException e){
            throw new GATKException(e.getMessage());
        }
    }

    /**
     * Converts a long string (typically from Spark) containing one FASTQ file contents, with line break "\n",
     * into a list of {@link FastqRecord}'s.
     */
    @VisibleForTesting
    public static List<FastqRecord> extractFASTQContents(final String fastqFileContents){
        final String[] contents = fastqFileContents.split("\n");
        final List<FastqRecord> result = new ArrayList<>(contents.length/4);
        for(int i=0; i<contents.length; i+=4){
            result.add(new FastqRecord(contents[i], contents[i+1], contents[i+2], contents[i+3]));
        }
        return result;
    }

    /**
     * Converts a list of FASTQRecords to a corresponding list of GATKRead's.
     * The input and output are guaranteed to have the same length.
     */
    public static List<GATKRead> convertToReads(final Iterable<FastqRecord> fastqRecords) {
        return StreamSupport.stream(fastqRecords.spliterator(), false).map(SVFastqUtils::convertToRead).collect(Collectors.toList());
    }

    /**
     * TODO: need to think more about this: do we want FastqRecord,
     * TODO:   or some custom record type that resembles the FASTQ, but with the original mapping and alignment information kept.
     * TODO:   in this case we need a parser. This is the parser, or the record type itself could provide such conversion to a GATKRead.
     * @param customFastq assumes the fastq record's read name contains two fields: the normal read name and "mapping=..." field, separated by a single space
     * @return reconstructed read will have an attribute "RC", and have value 1.
     */
    @VisibleForTesting
    static GATKRead convertToRead(final FastqRecord customFastq){
        final SAMRecord samRec = new SAMRecord(new SAMFileHeader());
        samRec.setReadString(customFastq.getReadString());
        samRec.setBaseQualities(SAMUtils.fastqToPhred(customFastq.getBaseQualityString()));
        // currently the FASTQ records template name line contains two fields separated by a space
        // see readToFastqRecord() in this file
        final String[] s = customFastq.getReadHeader().split(" ");
        samRec.setReadName(s[0].replace("@", ""));
        final String mapping = s[1].split("=")[1]; // the split would lead to a len 2 array {"mapping", chr:start;cigar} or {"mapping", "unmapped"}

        if (mapping.equals("unmapped")) {
            samRec.setReadUnmappedFlag(false);
        } else {
            final String[] t = mapping.split(":");
            samRec.setReferenceName(t[0]);
            final String[] loc = t[1].split(";");
            samRec.setAlignmentStart(Integer.valueOf(loc[0]));
            samRec.setCigarString(loc[1]);
        }

        final GATKRead read = SAMRecordToGATKReadAdapter.headerlessReadAdapter(samRec);
        read.setAttribute("RC", 1);
        return read;
    }
}
