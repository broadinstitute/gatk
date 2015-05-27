package org.broadinstitute.hellbender.utils.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.runners.BlockingDataflowPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.DirectPipelineRunner;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.util.GcsUtil;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.google.cloud.genomics.dataflow.utils.GenomicsDatasetOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.nio.channels.Channels;
import java.security.GeneralSecurityException;
import java.util.List;

/**
 * Utilities for working with google Dataflow
 *
 * Provides a a number of useful PTransforms and DoFns
 */
public final class DataflowUtils {

    public enum SaveDestination {
        LOCAL_DISK,
        CLOUD
    };

    private DataflowUtils(){} //prevent instantiation

    /**
     * a transform which will convert the input PCollection<I> to a PCollection<String> by calling toString() on each element
     * @return a Transform from I -> String
     */
    public static <I> PTransform<PCollection<? extends I>,PCollection<String>> convertToString(){
        return ParDo.of(
                new DoFn<I, String>() {
                  @Override
                  public void processElement(ProcessContext c) {
                      c.output(c.element().toString());
                  }
              });
    }

    /**
     * ingest local bam files from the file system and loads them into a PCollection<Read>
     * @param pipeline a configured Pipeline
     * @param intervals intervals to select reads from
     * @param bams paths to bam files to read from
     * @return a PCollection<Read> with all the reads the overlap the given intervals in the bams
     */
    public static PCollection<Read> getReadsFromLocalBams(final Pipeline pipeline, final List<SimpleInterval> intervals, final List<File> bams) {
        return pipeline.apply(Create.of(bams))
                .apply(ParDo.of(new LoadReadsFromFileFn(intervals)));
    }


    /**
     * get a transform that throws a specified exception
     */
    public static <I,O> PTransform<PCollection<? extends I>,PCollection<O>> throwException(Exception e){
        return ParDo.of(new ThrowExceptionFn<I, O>(e));
    }

    /**
     * throw a specified exception on execution
     */
    public static class ThrowExceptionFn<I,O> extends DoFn<I,O> {
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
    public static class LoadReadsFromFileFn extends DoFn<File, Read> {
        private final List<SimpleInterval> intervals;

        public LoadReadsFromFileFn(List<SimpleInterval> intervals) {
            this.intervals = intervals;
        }

        @Override
        public void processElement(ProcessContext c) {
            ReadsDataSource sams = new ReadsDataSource(c.element());
            sams.setIntervalsForTraversal(intervals);
            for (SAMRecord sam : sams) {
                c.output(ReadConverter.makeRead(sam));
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
     * @param fname the name of the destination, starting with "gs://" to save to GCS.
     * @returns SaveDestination.CLOUD if saved to GCS, SaveDestination.LOCAL_DISK otherwise.
     */
    public static <T> SaveDestination serializeSingleObject(PCollection<T> collection, String fname) {
        if (BucketUtils.isCloudStorageUrl(fname)) {
            saveSingleResultToGCS(collection, fname);
            return SaveDestination.CLOUD;
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
     * @param gcsDestPath the name of the destination (must start with "gs://").
     */
    public static <T> void saveSingleResultToGCS(final PCollection<T> collection, String gcsDestPath) {
        collection.apply(ParDo.named("save to " + gcsDestPath)
                .of(new DoFn<T, Void>() {
                    @Override
                    public void processElement(ProcessContext c) throws IOException, GeneralSecurityException {
                        GcsPath dest = GcsPath.fromUri(gcsDestPath);
                        GcsUtil gcsUtil = new GcsUtil.GcsUtilFactory().create(c.getPipelineOptions());
                        try (ObjectOutputStream out = new ObjectOutputStream(Channels.newOutputStream(gcsUtil.create(dest, "application/octet-stream")))) {
                            out.writeObject(c.element());
                        }
                    }
                }));
    }

}
