package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.util.IOUtil;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.Job;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.parquet.avro.AvroParquetOutputFormat;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.bdgenomics.adam.models.RecordGroupDictionary;
import org.bdgenomics.adam.models.SequenceDictionary;
import org.bdgenomics.formats.avro.AlignmentRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GATKReadToBDGAlignmentRecordConverter;
import org.broadinstitute.hellbender.utils.read.HeaderlessSAMRecordCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import org.disq_bio.disq.*;
import scala.Tuple2;

import java.io.IOException;
import java.util.Comparator;

/**
 * ReadsSparkSink writes GATKReads to a file. This code lifts from the HadoopGenomics/Hadoop-BAM
 * read writing code as well as from bigdatagenomics/adam.
 */
public final class ReadsSparkSink {

    private final static Logger logger = LogManager.getLogger(ReadsSparkSink.class);

    /**
     * writeReads writes rddReads to outputFile with header as the file header.
     * @param ctx the JavaSparkContext to write.
     * @param outputFile path to the output bam.
     * @param referenceFile path to the reference. required for cram output, otherwise may be null.
     * @param reads reads to write.
     * @param header the header to put at the top of the files
     * @param format should the output be a single file, sharded, ADAM, etc.
     */
    public static void writeReads(
            final JavaSparkContext ctx, final String outputFile, final String referenceFile, final JavaRDD<GATKRead> reads,
            final SAMFileHeader header, ReadsWriteFormat format) throws IOException {
        writeReads(ctx, outputFile, referenceFile, reads, header, format, 0, null);
    }

    /**
     * writeReads writes rddReads to outputFile with header as the file header.
     * @param ctx the JavaSparkContext to write.
     * @param outputFile path to the output bam.
     * @param referenceFile path to the reference. required for cram output, otherwise may be null.
     * @param reads reads to write.
     * @param header the header to put at the top of the files
     * @param format should the output be a single file, sharded, ADAM, etc.
     * @param numReducers the number of reducers to use when writing a single file. A value of zero indicates that the default
     *                    should be used.
     * @param outputPartsDir directory for temporary files for SINGLE output format, should be null for default value of filename + .output
     */
    public static void writeReads(
            final JavaSparkContext ctx, final String outputFile, final String referenceFile, final JavaRDD<GATKRead> reads,
            final SAMFileHeader header, ReadsWriteFormat format, final int numReducers, final String outputPartsDir) throws IOException {

        String absoluteOutputFile = BucketUtils.makeFilePathAbsolute(outputFile);
        String absoluteReferenceFile = referenceFile != null ?
                                        BucketUtils.makeFilePathAbsolute(referenceFile) :
                                        referenceFile;
      ReadsSparkSource.checkCramReference(ctx, absoluteOutputFile, absoluteReferenceFile);

        // The underlying reads are required to be in SAMRecord format in order to be
        // written out, so we convert them to SAMRecord explicitly here. If they're already
        // SAMRecords, this will effectively be a no-op. The SAMRecords will be headerless
        // for efficient serialization.
        final JavaRDD<SAMRecord> samReads = reads.map(read -> read.convertToSAMRecord(null));

        if (format == ReadsWriteFormat.SINGLE) {
            FileCardinalityWriteOption fileCardinalityWriteOption = FileCardinalityWriteOption.SINGLE;
            final String outputPartsDirectory = (outputPartsDir == null)? getDefaultPartsDirectory(outputFile)  : outputPartsDir;
            TempPartsDirectoryWriteOption tempPartsDirectoryWriteOption = new TempPartsDirectoryWriteOption(outputPartsDirectory);
            if (absoluteOutputFile.endsWith(BamFileIoUtils.BAM_FILE_EXTENSION) ||
                    absoluteOutputFile.endsWith(CramIO.CRAM_FILE_EXTENSION) ||
                    absoluteOutputFile.endsWith(IOUtil.SAM_FILE_EXTENSION)) {
                // don't specify a write option for format since it is inferred from the extension in the path
                writeReads(ctx, absoluteOutputFile, absoluteReferenceFile, samReads, header, numReducers, fileCardinalityWriteOption, tempPartsDirectoryWriteOption);
            } else {
                // default to BAM
                ReadsFormatWriteOption formatWriteOption = ReadsFormatWriteOption.BAM;
                writeReads(ctx, absoluteOutputFile, absoluteReferenceFile, samReads, header, numReducers, formatWriteOption, fileCardinalityWriteOption, tempPartsDirectoryWriteOption);
            }
        } else if (format == ReadsWriteFormat.SHARDED) {
            if (outputPartsDir!=null) {
                throw new  GATKException(String.format("You specified the bam output parts directory %s, but requested a sharded output format which does not use this option",outputPartsDir));
            }
            ReadsFormatWriteOption formatWriteOption = ReadsFormatWriteOption.BAM; // use BAM if output file is a directory
            FileCardinalityWriteOption fileCardinalityWriteOption = FileCardinalityWriteOption.MULTIPLE;
            writeReads(ctx, absoluteOutputFile, absoluteReferenceFile, samReads, header, numReducers, formatWriteOption, fileCardinalityWriteOption);
        } else if (format == ReadsWriteFormat.ADAM) {
            if (outputPartsDir!=null) {
                throw new  GATKException(String.format("You specified the bam output parts directory %s, but requested an ADAM output format which does not use this option",outputPartsDir));
            }
            writeReadsADAM(ctx, absoluteOutputFile, samReads, header);
        }
    }

    private static void writeReads(
            final JavaSparkContext ctx, final String outputFile, final String referenceFile, final JavaRDD<SAMRecord> reads,
            final SAMFileHeader header, final int numReducers, final WriteOption... writeOptions) throws IOException {

        final JavaRDD<SAMRecord> sortedReads = sortSamRecordsToMatchHeader(reads, header, numReducers);
        Broadcast<SAMFileHeader> headerBroadcast = ctx.broadcast(header);
        final JavaRDD<SAMRecord> sortedReadsWithHeader = sortedReads.map(read -> {
            read.setHeaderStrict(headerBroadcast.getValue());
            return read;
        });
        HtsjdkReadsRdd htsjdkReadsRdd = new HtsjdkReadsRdd(header, sortedReadsWithHeader);
        HtsjdkReadsRddStorage.makeDefault(ctx)
                .referenceSourcePath(referenceFile)
                .write(htsjdkReadsRdd, outputFile, writeOptions);
    }

    private static void writeReadsADAM(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<SAMRecord> reads,
            final SAMFileHeader header) throws IOException {
        final SequenceDictionary seqDict = SequenceDictionary.fromSAMSequenceDictionary(header.getSequenceDictionary());
        final RecordGroupDictionary readGroups = RecordGroupDictionary.fromSAMHeader(header);
        final JavaPairRDD<Void, AlignmentRecord> rddAlignmentRecords =
                reads.map(read -> {
                    read.setHeaderStrict(header);
                    AlignmentRecord alignmentRecord = GATKReadToBDGAlignmentRecordConverter.convert(read, seqDict, readGroups);
                    read.setHeaderStrict(null); // Restore the header to its previous state so as not to surprise the caller
                    return alignmentRecord;
                }).mapToPair(alignmentRecord -> new Tuple2<>(null, alignmentRecord));
        // instantiating a Job is necessary here in order to set the Hadoop Configuration...
        final Job job = Job.getInstance(ctx.hadoopConfiguration());
        // ...here, which sets a config property that the AvroParquetOutputFormat needs when writing data. Specifically,
        // we are writing the Avro schema to the Configuration as a JSON string. The AvroParquetOutputFormat class knows
        // how to translate objects in the Avro data model to the Parquet primitives that get written.
        AvroParquetOutputFormat.setSchema(job, AlignmentRecord.getClassSchema());
        deleteHadoopFile(outputFile, ctx.hadoopConfiguration());
        rddAlignmentRecords.saveAsNewAPIHadoopFile(
                outputFile, Void.class, AlignmentRecord.class, AvroParquetOutputFormat.class, job.getConfiguration());
    }

    private static void deleteHadoopFile(String fileToObliterate, Configuration conf) throws IOException {
        final Path pathToDelete = new Path(fileToObliterate);
        pathToDelete.getFileSystem(conf).delete(pathToDelete, true);
    }

    /**
     * Gets the default parts directory for a given file by appending .parts/ to the end of it
     */
    public static String getDefaultPartsDirectory(String file) {
        return file + ".parts/";
    }

    /**
     * Sorts the given reads according to the sort order in the header.
     * @param reads the reads to sort
     * @param header the header specifying the sort order,
     *               if the header specifies {@link SAMFileHeader.SortOrder#unsorted} or {@link SAMFileHeader.SortOrder#unknown}
     *               then no sort will be performed
     * @param numReducers the number of reducers to use; a value of 0 means use the default number of reducers
     * @return a sorted RDD of reads
     */
    private static JavaRDD<SAMRecord> sortSamRecordsToMatchHeader(final JavaRDD<SAMRecord> reads, final SAMFileHeader header, final int numReducers) {
        final Comparator<SAMRecord> comparator = getSAMRecordComparator(header);
        if ( comparator == null ) {
            return reads;
        } else {
            return SparkUtils.sortUsingElementsAsKeys(reads, comparator, numReducers);
        }
    }

    //Returns the comparator to use or null if no sorting is required.
    private static Comparator<SAMRecord> getSAMRecordComparator(final SAMFileHeader header) {
        switch (header.getSortOrder()){
            case coordinate: return new HeaderlessSAMRecordCoordinateComparator(header);
            //duplicate isn't supported because it doesn't work right on headerless SAMRecords
            case duplicate: throw new UserException.UnimplementedFeature("The sort order \"duplicate\" is not supported in Spark.");
            case queryname:
            case unsorted:   return header.getSortOrder().getComparatorInstance();
            default:         return null; //NOTE: javac warns if you have this (useless) default BUT it errors out if you remove this default.
        }
    }

}
