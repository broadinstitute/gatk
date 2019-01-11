package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.Job;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.parquet.avro.AvroParquetInputFormat;
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
import org.disq_bio.disq.HtsjdkReadsRdd;
import org.disq_bio.disq.HtsjdkReadsRddStorage;
import org.disq_bio.disq.HtsjdkReadsTraversalParameters;

import java.io.IOException;
import java.io.Serializable;
import java.util.List;
import java.util.Objects;

/** Loads the reads from disk either serially (using samReaderFactory) or in parallel using Hadoop-BAM.
 * The parallel code is a modified version of the example writing code from Hadoop-BAM.
 */
public final class ReadsSparkSource implements Serializable {
    private static final long serialVersionUID = 1L;

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
        try {
            String cramReferencePath = checkCramReference(ctx, readFileName, referencePath);
            HtsjdkReadsTraversalParameters<SimpleInterval> tp = traversalParameters == null ? null :
                    new HtsjdkReadsTraversalParameters<>(traversalParameters.getIntervalsForTraversal(), traversalParameters.traverseUnmappedReads());
            HtsjdkReadsRdd htsjdkReadsRdd = HtsjdkReadsRddStorage.makeDefault(ctx)
                    .splitSize((int) splitSize)
                    .validationStringency(validationStringency)
                    .referenceSourcePath(cramReferencePath)
                    .read(readFileName, tp);
            JavaRDD<GATKRead> reads = htsjdkReadsRdd.getReads()
                    .map(read -> (GATKRead) SAMRecordToGATKReadAdapter.headerlessReadAdapter(read))
                    .filter(Objects::nonNull);
            return fixPartitionsIfQueryGrouped(ctx, htsjdkReadsRdd.getHeader(), reads);
        } catch (IOException | IllegalArgumentException e) {
            throw new UserException("Failed to load reads from " + readFileName + "\n Caused by:" + e.getMessage(), e);
        }
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
            String cramReferencePath = checkCramReference(ctx, filePath, referencePath);
            return HtsjdkReadsRddStorage.makeDefault(ctx)
                    .validationStringency(validationStringency)
                    .referenceSourcePath(cramReferencePath)
                    .read(filePath)
                    .getHeader();
        } catch (IOException | IllegalArgumentException e) {
            throw new UserException("Failed to read bam header from " + filePath + "\n Caused by:" + e.getMessage(), e);
        }
    }

    /**
     * Check that for CRAM the reference is set to a file that exists and is not 2bit.
     * @return the <code>referencePath</code> or <code>null</code> if not CRAM
     */
    static String checkCramReference(final JavaSparkContext ctx, final String filePath, final String referencePath) {
        if (IOUtils.isCramFileName(filePath)) {
            if (referencePath == null) {
                throw new UserException.MissingReference("A reference is required for CRAM input");
            } else if (ReferenceTwoBitSparkSource.isTwoBit(referencePath)) { // htsjdk can't handle 2bit reference files
                throw new UserException("A 2bit file cannot be used as a CRAM file reference");
            } else {
                final Path refPath = new Path(referencePath);
                if (!SparkUtils.pathExists(ctx, refPath)) {
                    throw new UserException.MissingReference("The specified fasta file (" + referencePath + ") does not exist.");
                }
            }
            return referencePath;
        }
        return null;
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
