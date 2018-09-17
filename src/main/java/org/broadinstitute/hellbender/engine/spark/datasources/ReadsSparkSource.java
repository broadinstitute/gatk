package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.PathFilter;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.parquet.avro.AvroParquetInputFormat;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.bdgenomics.formats.avro.AlignmentRecord;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.TraversalParameters;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import org.seqdoop.hadoop_bam.*;
import org.seqdoop.hadoop_bam.util.SAMHeaderReader;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;

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
     * @param traversalParameters parameters controlling which reads to include. If <code>null</code> then all the reads (both mapped and unmapped) will be returned.
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String readFileName, final String referencePath, final TraversalParameters traversalParameters) {
        return getParallelReads(readFileName, referencePath, traversalParameters, 0);
    }


    /**
     * this is a hack to work around https://github.com/HadoopGenomics/Hadoop-BAM/issues/199
     *
     * fix the problem by explicitly sorting the input file splits
     */
    public static class SplitSortingSamInputFormat extends AnySAMInputFormat{
        @SuppressWarnings("unchecked")
        @Override
        public List<InputSplit> getSplits(JobContext job) throws IOException {
            final List<InputSplit> splits = super.getSplits(job);

            if( splits.stream().allMatch(split -> split instanceof FileVirtualSplit || split instanceof FileSplit)) {
                splits.sort(Comparator.comparing(split -> {
                    if (split instanceof FileVirtualSplit) {
                        return ((FileVirtualSplit) split).getPath().getName();
                    } else {
                        return ((FileSplit) split).getPath().getName();
                    }
                }));
            }

            return splits;
        }
    }


    /**
     * Loads Reads using Hadoop-BAM. For local files, bam must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param readFileName file to load
     * @param referencePath Reference path or null if not available. Reference is required for CRAM files.
     * @param traversalParameters parameters controlling which reads to include. If <code>null</code> then all the reads (both mapped and unmapped) will be returned.
     * @param splitSize maximum bytes of bam file to read into a single partition, increasing this will result in fewer partitions. A value of zero means
     *                  use the default split size (determined by the Hadoop input format, typically the size of one HDFS block).
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String readFileName, final String referencePath, final TraversalParameters traversalParameters, final long splitSize) {
        SAMFileHeader header = getHeader(readFileName, referencePath);

        // use the Hadoop configuration attached to the Spark context to maintain cumulative settings
        final Configuration conf = ctx.hadoopConfiguration();
        if (splitSize > 0) {
            conf.set("mapreduce.input.fileinputformat.split.maxsize", Long.toString(splitSize));
        }

        final JavaPairRDD<LongWritable, SAMRecordWritable> rdd2;

        setHadoopBAMConfigurationProperties(readFileName, referencePath);

        boolean isBam = IOUtils.isBamFileName(readFileName);
        if (isBam) {
            if (traversalParameters == null) {
                BAMInputFormat.unsetTraversalParameters(conf);
            } else {
                BAMInputFormat.setTraversalParameters(conf, traversalParameters.getIntervalsForTraversal(), traversalParameters.traverseUnmappedReads());
            }
        }

        rdd2 = ctx.newAPIHadoopFile(
                readFileName, SplitSortingSamInputFormat.class, LongWritable.class, SAMRecordWritable.class,
                conf);

        JavaRDD<GATKRead> reads= rdd2.map(v1 -> {
            SAMRecord sam = v1._2().get();
            if (isBam || samRecordOverlaps(sam, traversalParameters)) { // don't check overlaps for BAM since it is done by input format
                return (GATKRead) SAMRecordToGATKReadAdapter.headerlessReadAdapter(sam);
            }
            return null;
        }).filter(Objects::nonNull);

        return fixPartitionsIfQueryGrouped(ctx, header, reads);
    }

    private static JavaRDD<GATKRead> fixPartitionsIfQueryGrouped(JavaSparkContext ctx, SAMFileHeader header, JavaRDD<GATKRead> reads) {
        if( ReadUtils.isReadNameGroupedBam(header)) {
            return SparkUtils.putReadsWithTheSameNameInTheSamePartition(header, reads, ctx);
        } else {
            return reads;
        }
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, readFileName must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param readFileName file to load
     * @param referencePath Reference path or null if not available. Reference is required for CRAM files.
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String readFileName, final String referencePath) {
        return getParallelReads(readFileName, referencePath, 0);
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, readFileName must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param readFileName file to load
     * @param referencePath Reference path or null if not available. Reference is required for CRAM files.
     * @param splitSize maximum bytes of bam file to read into a single partition, increasing this will result in fewer partitions. A value of zero means
     *                  use the default split size (determined by the Hadoop input format, typically the size of one HDFS block).
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String readFileName, final String referencePath, int splitSize) {
        return getParallelReads(readFileName, referencePath, null /* all reads */, splitSize);
    }

    /**
     * Loads ADAM reads stored as Parquet.
     * @param inputPath path to the Parquet data
     * @return RDD of (ADAM-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getADAMReads(final String inputPath, final TraversalParameters traversalParameters, final SAMFileHeader header) throws IOException {
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
        JavaRDD<GATKRead> filteredRdd = readsRdd.filter(record -> samRecordOverlaps(record.convertToSAMRecord(header), traversalParameters));

        return fixPartitionsIfQueryGrouped(ctx, header, filteredRdd);
    }

    /**
     * Loads the header using Hadoop-BAM.
     * @param filePath path to the bam.
     * @param referencePath Reference path or null if not available. Reference is required for CRAM files.
     * @return the header for the bam.
     */
    public SAMFileHeader getHeader(final String filePath, final String referencePath) {
        // GCS case
        if (BucketUtils.isCloudStorageUrl(filePath)) {
            try (ReadsDataSource readsDataSource = new ReadsDataSource(IOUtils.getPath(filePath))) {
                return readsDataSource.getHeader();
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
        } catch (IOException | IllegalArgumentException e) {
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
                if ( ReferenceTwoBitSparkSource.isTwoBit(referenceName)) { // htsjdk can't handle 2bit reference files
                    throw new UserException("A 2bit file cannot be used as a CRAM file reference");
                }
                else { // Hadoop-BAM requires the reference to be a URI, including scheme
                    final Path refPath = new Path(referenceName);
                    if (!SparkUtils.pathExists(ctx, refPath)) {
                        throw new UserException.MissingReference("The specified fasta file (" + referenceName + ") does not exist.");
                    }
                    final String referenceURI = null == refPath.toUri().getScheme() ?
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
    private static boolean samRecordOverlaps(final SAMRecord record, final TraversalParameters traversalParameters ) {
        if (traversalParameters == null) {
            return true;
        }
        if (traversalParameters.traverseUnmappedReads() && record.getReadUnmappedFlag() && record.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
            return true; // include record if unmapped records should be traversed and record is unmapped
        }
        List<SimpleInterval> intervals = traversalParameters.getIntervalsForTraversal();
        if (intervals == null || intervals.isEmpty()) {
            return false; // no intervals means 'no mapped reads'
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
