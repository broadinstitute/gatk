package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.api.services.storage.Storage;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import htsjdk.samtools.*;
import htsjdk.samtools.cram.build.CramIO;
import org.apache.commons.lang3.StringUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.PathFilter;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.Job;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.parquet.avro.AvroParquetInputFormat;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.bdgenomics.formats.avro.AlignmentRecord;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.BDGAlignmentRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadConstants;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.seqdoop.hadoop_bam.AnySAMInputFormat;
import org.seqdoop.hadoop_bam.BAMInputFormat;
import org.seqdoop.hadoop_bam.CRAMInputFormat;
import org.seqdoop.hadoop_bam.SAMRecordWritable;
import org.seqdoop.hadoop_bam.util.SAMHeaderReader;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.List;

/** Loads the reads from disk either serially (using samReaderFactory) or in parallel using Hadoop-BAM.
 * The parallel code is a modified version of the example writing code from Hadoop-BAM.
 */
public final class ReadsSparkSource implements Serializable {
    private static final long serialVersionUID = 1L;
    private static final String HADOOP_PART_PREFIX = "part-";

    private transient final JavaSparkContext ctx;
    private ValidationStringency validationStringency = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;

    private static final Logger logger = LogManager.getLogger(ReadsSparkSource.class);

    public ReadsSparkSource(final JavaSparkContext ctx) { this.ctx = ctx; }

    public ReadsSparkSource(final JavaSparkContext ctx, final ValidationStringency validationStringency)
    {
        this.ctx = ctx;
        this.validationStringency = validationStringency;
    }


    /**
     * Loads Reads using Hadoop-BAM. For local files, readFileName must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param readFileName file to load
     * @param referencePath Reference path or null if not available. Reference is required for CRAM files.
     * @param intervals intervals of reads to include.
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String readFileName, final String referencePath, final List<SimpleInterval> intervals) {
        return getParallelReads(readFileName, referencePath, intervals, 0);
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, bam must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param readFileName file to load
     * @param referencePath Reference path or null if not available. Reference is required for CRAM files.
     * @param intervals intervals of reads to include.
     * @param splitSize maximum bytes of bam file to read into a single partition, increasing this will result in fewer partitions. A value of zero means
     *                  use the default split size (determined by the Hadoop input format, typically the size of one HDFS block).
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String readFileName, final String referencePath, final List<SimpleInterval> intervals, final long splitSize) {
        // use the Hadoop configuration attached to the Spark context to maintain cumulative settings
        final Configuration conf = ctx.hadoopConfiguration();
        if (splitSize > 0) {
            conf.set("mapreduce.input.fileinputformat.split.maxsize", Long.toString(splitSize));
        }

        final JavaPairRDD<LongWritable, SAMRecordWritable> rdd2;

        setHadoopBAMConfigurationProperties(readFileName, referencePath);

        boolean isBam = IOUtils.isBamFileName(readFileName);
        if (isBam && intervals != null && !intervals.isEmpty()) {
            BAMInputFormat.setIntervals(conf, intervals);
        } else {
            conf.unset(BAMInputFormat.INTERVALS_PROPERTY);
        }

        rdd2 = ctx.newAPIHadoopFile(
                    readFileName, AnySAMInputFormat.class, LongWritable.class, SAMRecordWritable.class,
                    conf);

        return rdd2.map(v1 -> {
            SAMRecord sam = v1._2().get();
            if (isBam || samRecordOverlaps(sam, intervals)) { // don't check overlaps for BAM since it is done by input format
                return (GATKRead) SAMRecordToGATKReadAdapter.headerlessReadAdapter(sam);
            }
            return null;
        }).filter(v1 -> v1 != null);
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, readFileName must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam. This excludes unmapped reads.
     * @param readFileName file to load
     * @param referencePath Reference path or null if not available. Reference is required for CRAM files.
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String readFileName, final String referencePath) {
        return getParallelReads(readFileName, referencePath, 0);
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, readFileName must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam. This excludes unmapped reads.
     * @param readFileName file to load
     * @param referencePath Reference path or null if not available. Reference is required for CRAM files.
     * @param splitSize maximum bytes of bam file to read into a single partition, increasing this will result in fewer partitions. A value of zero means
     *                  use the default split size (determined by the Hadoop input format, typically the size of one HDFS block).
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String readFileName, final String referencePath, int splitSize) {
        final SAMFileHeader readsHeader = getHeader(readFileName, referencePath, null);
        List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());
        return getParallelReads(readFileName, referencePath, intervals, splitSize);
    }

    /**
     * Loads ADAM reads stored as Parquet.
     * @param inputPath path to the Parquet data
     * @return RDD of (ADAM-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getADAMReads(final String inputPath, final List<SimpleInterval> intervals, final SAMFileHeader header) throws IOException {
        Job job = Job.getInstance(ctx.hadoopConfiguration());
        AvroParquetInputFormat.setAvroReadSchema(job, AlignmentRecord.getClassSchema());
        Broadcast<SAMFileHeader> bHeader;
        if (header == null) {
            bHeader= ctx.broadcast(null);
        } else {
            bHeader = ctx.broadcast(header);
        }
        @SuppressWarnings("unchecked")
        JavaRDD<AlignmentRecord> recordsRdd = ctx.newAPIHadoopFile(
                inputPath, AvroParquetInputFormat.class, Void.class, AlignmentRecord.class, job.getConfiguration())
                .values();
        JavaRDD<GATKRead> readsRdd = recordsRdd.map(record -> new BDGAlignmentRecordToGATKReadAdapter(record, bHeader.getValue()));
        JavaRDD<GATKRead> filteredRdd = readsRdd.filter(record -> samRecordOverlaps(record.convertToSAMRecord(header), intervals));
        return filteredRdd;
    }

    /**
     * Loads the header using Hadoop-BAM.
     * @param filePath path to the bam.
     * @param referencePath Reference path or null if not available. Reference is required for CRAM files.
     * @param auth authentication information if using GCS.
     * @return the header for the bam.
     */
    public SAMFileHeader getHeader(final String filePath, final String referencePath, final AuthHolder auth) {
        // GCS case
        if (BucketUtils.isCloudStorageUrl(filePath)) {
            try {
                Storage.Objects storageClient = auth.makeStorageClient();
                try (final SamReader reader = BAMIO.openBAM(storageClient, filePath, validationStringency)) {
                    return reader.getFileHeader();
                }
            } catch (Exception e) {
                throw new UserException("Failed to read bam header from " + filePath + "\n Caused by:" + e.getMessage(), e);
            }
        }

        // local file or HDFs case
        try {
            Path path = new Path(filePath);
            FileSystem fs = path.getFileSystem(ctx.hadoopConfiguration());
            if (fs.isDirectory(path)) {
                FileStatus[] bamFiles = fs.listStatus(path, new PathFilter() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public boolean accept(Path path) {
                        return path.getName().startsWith(HADOOP_PART_PREFIX);
                    }
                });
                if (bamFiles.length == 0) {
                    throw new UserException("No BAM files to load header from in: " + path);
                }
                path = bamFiles[0].getPath(); // Hadoop-BAM writes the same header to each shard, so use the first one
            }
            setHadoopBAMConfigurationProperties(filePath, referencePath);
            return SAMHeaderReader.readSAMHeaderFrom(path, ctx.hadoopConfiguration());
        } catch (IOException e) {
            throw new UserException("Failed to read bam header from " + filePath + "\n Caused by:" + e.getMessage(), e);
        }
    }

    /**
     * Propagate any values that need to be passed to Hadoop-BAM through configuration properties:
     *
     *   - the validation stringency property is always set using the current value of the
     *     validationStringency field
     *   - if the input file is a CRAM file, the reference value will also be set, and must be a URI
     *     which includes a scheme. if no scheme is provided a "file://" scheme will be used. for
     *     non-CRAM input the reference may be null.
     *   - if the input file is not CRAM, the reference property is *unset* to prevent Hadoop-BAM
     *     from passing a stale value through to htsjdk when multiple read calls are made serially
     *     with different inputs but the same Spark context
     */
    private void setHadoopBAMConfigurationProperties(final String inputName, final String referenceName) {
        // use the Hadoop configuration attached to the Spark context to maintain cumulative settings
        final Configuration conf = ctx.hadoopConfiguration();
        conf.set(SAMHeaderReader.VALIDATION_STRINGENCY_PROPERTY, validationStringency.name());

        if (!IOUtils.isCramFileName(inputName)) {
            // only set the reference for CRAM input
            conf.unset(CRAMInputFormat.REFERENCE_SOURCE_PATH_PROPERTY);
        }
        else {
            if (null == referenceName) {
                throw new UserException.MissingReference("A reference is required for CRAM input");
            }
            else {
                if (ReferenceTwoBitSource.isTwoBit(referenceName)) { // htsjdk can't handle 2bit reference files
                    throw new UserException("A 2bit file cannot be used as a CRAM file reference");
                }
                else { // Hadoop-BAM requires the reference to be a URI, including scheme
                    String referenceURI = null ==  new Path(referenceName).toUri().getScheme() ?
                            "file://" + new File(referenceName).getAbsolutePath() :
                            referenceName;
                    conf.set(CRAMInputFormat.REFERENCE_SOURCE_PATH_PROPERTY, referenceURI);
                }
            }
        }
    }

    /**
     * Tests if a given SAMRecord overlaps any interval in a collection. This is only used as a fallback option for
     * formats that don't support query-by-interval natively at the Hadoop-BAM layer.
     */
    //TODO: use IntervalsSkipList, see https://github.com/broadinstitute/gatk/issues/1531
    private static boolean samRecordOverlaps(final SAMRecord record, final List<SimpleInterval> intervals ) {
        if (intervals == null || intervals.isEmpty()) {
            return true;
        }
        for (SimpleInterval interval : intervals) {
            if (record.getReadUnmappedFlag() && record.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START) {
                // This follows the behavior of htsjdk's SamReader which states that "an unmapped read will be returned
                // by this call if it has a coordinate for the purpose of sorting that is in the query region".
                int start = record.getAlignmentStart();
                return interval.getStart() <= start && interval.getEnd() >= start;
            } else  if (interval.overlaps(record)) {
                return true;
            }
        }
        return false;
    }
}
