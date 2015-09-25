package org.broadinstitute.hellbender.utils.dataflow;

import com.cloudera.dataflow.hadoop.HadoopIO;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.hadoop.io.LongWritable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.coders.ReadContextDataCoder;
import org.broadinstitute.hellbender.engine.dataflow.coders.UUIDCoder;
import org.broadinstitute.hellbender.engine.dataflow.coders.VariantCoder;
import org.broadinstitute.hellbender.engine.dataflow.datasources.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.MinimalVariant;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.seqdoop.hadoop_bam.AnySAMInputFormat;
import org.seqdoop.hadoop_bam.SAMRecordWritable;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.List;
import java.util.UUID;

/**
 * Utilities for working with google Dataflow
 *
 * Provides a a number of useful PTransforms and DoFns
 */
public final class DataflowUtils {

    private static final Logger logger = LogManager.getLogger(DataflowUtils.class);

    public enum SaveDestination {
        LOCAL_DISK,
        CLOUD,
        HDFS
    }

    private DataflowUtils(){} //prevent instantiation

    /**
     * Standard method for registering all coders needed by the GATK
     *
     * @param p pipeline for which to register coders
     */
    public static void registerGATKCoders( final Pipeline p ) {
        DataflowWorkarounds.registerGenomicsCoders(p);
        p.getCoderRegistry().registerCoder(GATKRead.class, new GATKReadCoder());
        p.getCoderRegistry().registerCoder(GoogleGenomicsReadToGATKReadAdapter.class, GoogleGenomicsReadToGATKReadAdapter.CODER);
        p.getCoderRegistry().registerCoder(SAMRecordToGATKReadAdapter.class, SerializableCoder.of(SAMRecordToGATKReadAdapter.class));
        p.getCoderRegistry().registerCoder(SimpleInterval.class, SerializableCoder.of(SimpleInterval.class));
        p.getCoderRegistry().registerCoder(UUID.class, UUIDCoder.CODER);
        p.getCoderRegistry().registerCoder(Variant.class, new VariantCoder());
        p.getCoderRegistry().registerCoder(VariantContextVariantAdapter.class, SerializableCoder.of(VariantContextVariantAdapter.class));
        p.getCoderRegistry().registerCoder(MinimalVariant.class, SerializableCoder.of(MinimalVariant.class));
        p.getCoderRegistry().registerCoder(RefAPISource.class, SerializableCoder.of(RefAPISource.class));
        p.getCoderRegistry().registerCoder(ReferenceDataflowSource.class, SerializableCoder.of(ReferenceDataflowSource.class));
        p.getCoderRegistry().registerCoder(ReferenceBases.class, SerializableCoder.of(ReferenceBases.class));
        p.getCoderRegistry().registerCoder(ReadContextData.class, new ReadContextDataCoder());
        p.getCoderRegistry().registerCoder(ReferenceShard.class, ReferenceShard.CODER);
        p.getCoderRegistry().registerCoder(VariantShard.class, VariantShard.CODER);
    }

    /**
     * a transform which will convert the input PCollection<I> to a PCollection<String> by calling toString() on each element
     * @return a Transform from I -> String
     */
    public static <I> PTransform<PCollection<? extends I>,PCollection<String>> convertToString(){
        return ParDo.of(
                new DoFn<I, String>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement( ProcessContext c ) {
                        c.output(c.element().toString());
                    }
                });
    }

    /**
     * a transform that prints the contents to standard out. This is meant for local testing.
     */
    public static class PrintCollection<T> extends DoFn<T, Void> {
        private final String name;
        public PrintCollection(String name) {
            this.name = name;
        }

        private static final long serialVersionUID = 1L;
        @Override
        public void processElement(ProcessContext c) throws Exception {
            T i = c.element();
            System.out.println(name + ": " + i);
        }
    }

    /**
     * ingest local bam files from the file system and loads them into a PCollection<GATKRead>
     * @param pipeline a configured Pipeline
     * @param intervals intervals to select reads from
     * @param bams paths to bam files to read from
     * @return a PCollection<GATKRead> with all the reads the overlap the given intervals in the bams
     */
    public static PCollection<GATKRead> getReadsFromLocalBams(final Pipeline pipeline, final List<SimpleInterval> intervals, final List<File> bams) {
        return getReadsFromLocalBams(pipeline, intervals, ValidationStringency.SILENT, bams);
    }

    /**
     * ingest local bam files from the file system and loads them into a PCollection<Read>
     * @param pipeline a configured Pipeline
     * @param intervals intervals to select reads from
     * @param stringency stringency of the input validation checks
     * @param bams paths to bam files to read from
     * @return a PCollection<GATKRead> with all the reads the overlap the given intervals in the bams
     */
    public static PCollection<GATKRead> getReadsFromLocalBams(final Pipeline pipeline, final List<SimpleInterval> intervals, final ValidationStringency stringency, final List<File> bams) {
        return pipeline.apply(Create.of(bams))
                .apply(ParDo.of(new LoadReadsFromFileFn(intervals, stringency)));
    }

    /**
     * Ingest a BAM file from a Hadoop file system and loads into a
     * <code>PCollection<Read></code>.
     * @param pipeline a configured Pipeline
     * @param intervals intervals to select reads from
     * @param bam Hadoop file path to read from
     * @return a <code>PCollection<Read></code> with all the reads that overlap the
     * given intervals in the BAM file
     */
    @SuppressWarnings("unchecked")
    public static PCollection<GATKRead> getReadsFromHadoopBam(final Pipeline pipeline, final List<SimpleInterval> intervals, final ValidationStringency stringency, final String bam) {
        PCollection<KV<LongWritable, SAMRecordWritable>> input = pipeline.apply(
                HadoopIO.Read.from(bam, AnySAMInputFormat.class, LongWritable.class, SAMRecordWritable.class));
        return input.apply(ParDo.of(new DoFn<KV<LongWritable, SAMRecordWritable>, GATKRead>() {
            private static final long serialVersionUID = 1L;
            @Override
            public void processElement( ProcessContext c ) throws Exception {
                SAMRecord sam = c.element().getValue().get();
                if ( samRecordOverlaps(sam, intervals) ) {
                    try {
                        c.output(new SAMRecordToGATKReadAdapter(sam));
                    }
                    catch ( SAMException e ) {
                        if ( stringency == ValidationStringency.STRICT ) {
                            throw e;
                        }
                        else if ( stringency == ValidationStringency.LENIENT ) {
                            logger.info("getReadsFromHadoopBam: " + e.getMessage());
                        }
                        // do nothing if silent
                    }
                }
            }
        }));
    }

    /**
     * Tests if a given SAMRecord overlaps any interval in a collection.
     */
    //TODO: remove this method when https://github.com/broadinstitute/hellbender/issues/559 is fixed
    private static boolean samRecordOverlaps( SAMRecord record, List<SimpleInterval> intervals ) {
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

    /**
     * get a transform that throws a specified exception
     */
    public static <I,O> PTransform<PCollection<? extends I>,PCollection<O>> throwException(Exception e){
        return ParDo.of(new ThrowExceptionFn<>(e));
    }

    /**
     * throw a specified exception on execution
     */
    public static class ThrowExceptionFn<I,O> extends DoFn<I,O> {
        private static final long serialVersionUID = 1L;
        private final Exception e;

        public ThrowExceptionFn(Exception e){
            this.e = e;
        }

        @Override
        public void processElement(ProcessContext c) throws Exception {
            throw e;
        }
    }

    /**
     * Read a bam file and output each of the reads in it
     */
    public static class LoadReadsFromFileFn extends DoFn<File, GATKRead> {
        private static final long serialVersionUID = 1L;
        private static final Logger logger = LogManager.getLogger(LoadReadsFromFileFn.class);

        private final List<SimpleInterval> intervals;
        private final ValidationStringency stringency;

        public LoadReadsFromFileFn(List<SimpleInterval> intervals, final ValidationStringency stringency) {
            this.intervals = intervals;
            this.stringency = stringency;
        }

        @Override
        public void processElement(ProcessContext c) {
            final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(stringency);

            ReadsDataSource bam = new ReadsDataSource(c.element(),samReaderFactory);
            bam.setIntervalsForTraversal(intervals);
            for ( GATKRead read : bam ) {
                c.output(read);
            }
        }
    }


    /**
     * Serializes the collection's single object to the specified file.
     *
     * Of course if you run on the cloud and specify a local path, the file will be saved
     * on a cloud worker, which may not be very useful.
     *
     * @param collection A collection with a single serializable object to save.
     * @param fname the name of the destination, starting with "gs://" to save to GCS, or "hdfs://" to save to HDFS.
     * @returns SaveDestination.CLOUD if saved to GCS, SaveDestination.HDFS if saved to HDFS,
     * SaveDestination.LOCAL_DISK otherwise.
     */
    public static <T> SaveDestination serializeSingleObject(PCollection<T> collection, String fname) {
        if (BucketUtils.isCloudStorageUrl(fname)) {
            saveSingleResultToRemoteStorage(collection, fname);
            return SaveDestination.CLOUD;
        } else if (BucketUtils.isHadoopUrl(fname)) {
            saveSingleResultToRemoteStorage(collection, fname);
            return SaveDestination.HDFS;
        } else {
            saveSingleResultToLocalDisk(collection, fname);
            return SaveDestination.LOCAL_DISK;
        }
    }

    /**
     * Serializes the collection's single object to the specified file.
     *
     * @param collection A collection with a single serializable object to save.
     * @param fname the name of the destination.
     */
    public static <T> void saveSingleResultToLocalDisk(PCollection<T> collection, String fname) {
        collection.apply(ParDo
                .named("save to " + fname)
                .of(new DoFn<T, Void>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement(ProcessContext c) throws IOException {
                        T obj = c.element();
                        try (ObjectOutputStream os = new ObjectOutputStream(new FileOutputStream(fname))) {
                            os.writeObject(obj);
                        }
                    }
                }));
    }

    /**
     * Serializes the collection's single object to the specified file.
     *
     * @param collection A collection with a single serializable object to save.
     * @param destPath the name of the destination (must start with "gs://" or "hdfs://").
     */
    public static <T> void saveSingleResultToRemoteStorage(final PCollection<T> collection, String destPath) {
        collection.apply(ParDo.named("save to " + destPath)
                .of(new DoFn<T, Void>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement(ProcessContext c) throws IOException {
                        try (ObjectOutputStream out = new ObjectOutputStream(BucketUtils.createFile(destPath, c.getPipelineOptions()))) {
                            out.writeObject(c.element());
                        }
                    }
                }));
    }

}
