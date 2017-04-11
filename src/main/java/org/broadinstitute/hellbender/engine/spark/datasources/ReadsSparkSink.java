package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.cram.build.CramIO;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapred.FileAlreadyExistsException;
import org.apache.hadoop.mapreduce.*;
import org.apache.parquet.avro.AvroParquetOutputFormat;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.bdgenomics.adam.models.RecordGroupDictionary;
import org.bdgenomics.adam.models.SequenceDictionary;
import org.bdgenomics.formats.avro.AlignmentRecord;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GATKReadToBDGAlignmentRecordConverter;
import org.broadinstitute.hellbender.utils.read.HeaderlessSAMRecordCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import org.seqdoop.hadoop_bam.*;
import org.seqdoop.hadoop_bam.util.SAMFileMerger;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;

/**
 * ReadsSparkSink writes GATKReads to a file. This code lifts from the HadoopGenomics/Hadoop-BAM
 * read writing code as well as from bigdatagenomics/adam.
 */
public final class ReadsSparkSink {

    // Output format class for writing BAM files through saveAsNewAPIHadoopFile. Must be public.
    public static class SparkBAMOutputFormat extends KeyIgnoringBAMOutputFormat<NullWritable> {
        public static SAMFileHeader bamHeader = null;

        public static void setHeader(final SAMFileHeader header) {
            bamHeader = header;
        }

        @Override
        public RecordWriter<NullWritable, SAMRecordWritable> getRecordWriter(TaskAttemptContext ctx) throws IOException {
            setSAMHeader(bamHeader);
            // use BAM extension for part files since they are valid BAM files
            return getRecordWriter(ctx, getDefaultWorkFile(ctx, BamFileIoUtils.BAM_FILE_EXTENSION));
        }

        @Override
        public void checkOutputSpecs(JobContext job) throws IOException {
            try {
                super.checkOutputSpecs(job);
            } catch (FileAlreadyExistsException e) {
                // delete existing files before overwriting them
                final Path outDir = getOutputPath(job);
                outDir.getFileSystem(job.getConfiguration()).delete(outDir, true);
            }
        }
    }

    public static class SparkHeaderlessBAMOutputFormat extends SparkBAMOutputFormat {
        public SparkHeaderlessBAMOutputFormat() {
            setWriteHeader(false);
        }
    }

    // Output format class for writing CRAM files through saveAsNewAPIHadoopFile. Must be public.
    public static class SparkCRAMOutputFormat extends KeyIgnoringCRAMOutputFormat<NullWritable> {
        public static SAMFileHeader bamHeader = null;

        public static void setHeader(final SAMFileHeader header) {
            bamHeader = header;
        }

        @Override
        public RecordWriter<NullWritable, SAMRecordWritable> getRecordWriter(TaskAttemptContext ctx) throws IOException {
            setSAMHeader(bamHeader);
            // use CRAM extension for part files since they are valid CRAM files
            return getRecordWriter(ctx, getDefaultWorkFile(ctx, CramIO.CRAM_FILE_EXTENSION));
        }

        @Override
        public void checkOutputSpecs(JobContext job) throws IOException {
            try {
                super.checkOutputSpecs(job);
            } catch (FileAlreadyExistsException e) {
                // delete existing files before overwriting them
                Path outDir = getOutputPath(job);
                outDir.getFileSystem(job.getConfiguration()).delete(outDir, true);
            }
        }
    }

    public static class SparkHeaderlessCRAMOutputFormat extends SparkCRAMOutputFormat {
        public SparkHeaderlessCRAMOutputFormat() {
            setWriteHeader(false);
        }
    }

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
        writeReads(ctx, outputFile, referenceFile, reads, header, format, 0);
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
     */
    public static void writeReads(
            final JavaSparkContext ctx, final String outputFile, final String referenceFile, final JavaRDD<GATKRead> reads,
            final SAMFileHeader header, ReadsWriteFormat format, final int numReducers) throws IOException {

        SAMFormat samOutputFormat = IOUtils.isCramFileName(outputFile) ? SAMFormat.CRAM : SAMFormat.BAM;

        String absoluteOutputFile = BucketUtils.makeFilePathAbsolute(outputFile);
        String absoluteReferenceFile = referenceFile != null ?
                                        BucketUtils.makeFilePathAbsolute(referenceFile) :
                                        referenceFile;
        setHadoopBAMConfigurationProperties(ctx, absoluteOutputFile, absoluteReferenceFile, format);

        // The underlying reads are required to be in SAMRecord format in order to be
        // written out, so we convert them to SAMRecord explicitly here. If they're already
        // SAMRecords, this will effectively be a no-op. The SAMRecords will be headerless
        // for efficient serialization.
        final JavaRDD<SAMRecord> samReads = reads.map(read -> read.convertToSAMRecord(null));

        if (format == ReadsWriteFormat.SINGLE) {
            writeReadsSingle(ctx, absoluteOutputFile, absoluteReferenceFile, samOutputFormat, samReads, header, numReducers);
        } else if (format == ReadsWriteFormat.SHARDED) {
            saveAsShardedHadoopFiles(ctx, absoluteOutputFile, absoluteReferenceFile, samOutputFormat, samReads, header, true);
        } else if (format == ReadsWriteFormat.ADAM) {
            writeReadsADAM(ctx, absoluteOutputFile, samReads, header);
        }
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

    private static void saveAsShardedHadoopFiles(
            final JavaSparkContext ctx, final String outputFile, final String referenceFile,
            final SAMFormat samOutputFormat, final JavaRDD<SAMRecord> reads, final SAMFileHeader header,
            final boolean writeHeader) throws IOException {
        // Set the static header on the driver thread.
        if (samOutputFormat == SAMFormat.CRAM) {
            SparkCRAMOutputFormat.setHeader(header);
        } else {
            SparkBAMOutputFormat.setHeader(header);
        }

        final Broadcast<SAMFileHeader> headerBroadcast = ctx.broadcast(header);

        // SparkBAM/CRAMOutputFormat are static classes, so we need to copy the header to each worker then call
        final JavaRDD<SAMRecord> readsRDD = setHeaderForEachPartition(reads, samOutputFormat, headerBroadcast);

        // The expected format for writing is JavaPairRDD where the key is ignored and the value is SAMRecordWritable.
        final JavaPairRDD<SAMRecord, SAMRecordWritable> rddSamRecordWriteable = pairReadsWithSAMRecordWritables(headerBroadcast, readsRDD);

        rddSamRecordWriteable.saveAsNewAPIHadoopFile(outputFile, SAMRecord.class, SAMRecordWritable.class, getOutputFormat(samOutputFormat, writeHeader), ctx.hadoopConfiguration());
    }

    /**
     * SparkBAM/CRAMOutputFormat has a static header value which must be set on each executor.
     */
    private static JavaRDD<SAMRecord> setHeaderForEachPartition(final JavaRDD<SAMRecord> reads, final SAMFormat samOutputFormat, final Broadcast<SAMFileHeader> headerBroadcast) {
        if (samOutputFormat == SAMFormat.CRAM) {
            return reads.mapPartitions(readIterator -> {
                SparkCRAMOutputFormat.setHeader(headerBroadcast.getValue());
                return readIterator;
            });
        }
        else {
            return reads.mapPartitions(readIterator -> {
                SparkBAMOutputFormat.setHeader(headerBroadcast.getValue());
                return readIterator;
            });
        }
    }

    private static void writeReadsSingle(
            final JavaSparkContext ctx, final String outputFile, final String referenceFile, final SAMFormat samOutputFormat, final JavaRDD<SAMRecord> reads,
            final SAMFileHeader header, final int numReducers) throws IOException {

        final JavaRDD<SAMRecord> sortedReads = SparkUtils.sortReads(reads, header, numReducers);
        final String outputPartsDirectory = outputFile + ".parts";
        saveAsShardedHadoopFiles(ctx, outputPartsDirectory, referenceFile, samOutputFormat, sortedReads,  header, false);
        SAMFileMerger.mergeParts(outputPartsDirectory, outputFile, samOutputFormat, header);
    }

    private static Class<? extends OutputFormat<NullWritable, SAMRecordWritable>> getOutputFormat(final SAMFormat samFormat, final boolean writeHeader) {
        if (samFormat == SAMFormat.CRAM) {
            return writeHeader ? SparkCRAMOutputFormat.class : SparkHeaderlessCRAMOutputFormat.class;
        }
        else {
            return writeHeader ? SparkBAMOutputFormat.class : SparkHeaderlessBAMOutputFormat.class;
        }
    }

    private static JavaPairRDD<SAMRecord, SAMRecordWritable> pairReadsWithSAMRecordWritables(Broadcast<SAMFileHeader> headerBroadcast, JavaRDD<SAMRecord> records) {
        return records.mapToPair(read -> {
            read.setHeaderStrict(headerBroadcast.getValue());
            final SAMRecordWritable samRecordWritable = new SAMRecordWritable();
            samRecordWritable.set(read);
            return new Tuple2<>(read, samRecordWritable);
        });
    }

    private static void deleteHadoopFile(String fileToObliterate, Configuration conf) throws IOException {
        final Path pathToDelete = new Path(fileToObliterate);
        pathToDelete.getFileSystem(conf).delete(pathToDelete, true);
    }

    /**
     * Propagate any values that need to be passed to Hadoop-BAM through configuration properties:
     *
     *   - if the output file is a CRAM file, the reference value, which must be a URI which includes
     *     a scheme, will also be set
     *   - if the output file is not CRAM, the reference property is *unset* to prevent Hadoop-BAM
     *     from passing a stale value through to htsjdk when multiple calls are made serially
     *     with different outputs but the same Spark context
     */
    private static void setHadoopBAMConfigurationProperties(final JavaSparkContext ctx, final String outputName,
                                                            final String referenceName, final ReadsWriteFormat format) {
        final Configuration conf = ctx.hadoopConfiguration();

        if (!IOUtils.isCramFileName(outputName)) { // only set the reference for CRAM output
            conf.unset(CRAMInputFormat.REFERENCE_SOURCE_PATH_PROPERTY);
            if (format == ReadsWriteFormat.SINGLE && IOUtils.isBamFileName(outputName)) {
                conf.setBoolean(BAMOutputFormat.WRITE_SPLITTING_BAI, true);
            } else {
                conf.setBoolean(BAMOutputFormat.WRITE_SPLITTING_BAI, false);
            }
        }
        else {
            if (null == referenceName) {
                throw new UserException.MissingReference("A reference is required for CRAM output");
            }
            else {
                if (ReferenceTwoBitSource.isTwoBit(referenceName)) { // htsjdk can't handle 2bit reference files
                    throw new UserException("A 2bit file cannot be used as a CRAM file reference");
                }
                else { // Hadoop-BAM requires the reference to be a URI, including scheme
                    String referenceURI =
                            null == new Path(referenceName).toUri().getScheme() ?
                                "file://" + new File(referenceName).getAbsolutePath() :
                                referenceName;
                    conf.set(CRAMInputFormat.REFERENCE_SOURCE_PATH_PROPERTY, referenceURI);
                }
            }
        }
    }

}
