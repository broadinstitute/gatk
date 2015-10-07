package org.broadinstitute.hellbender.utils.dataflow;


import com.cloudera.dataflow.hadoop.HadoopIO;
import com.cloudera.dataflow.hadoop.WritableCoder;
import com.cloudera.dataflow.spark.ShardNameTemplateAware;
import com.cloudera.dataflow.spark.ShardNameTemplateHelper;
import com.cloudera.dataflow.spark.SparkPipelineRunner;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.*;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.util.SerializableUtils;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import org.apache.hadoop.conf.Configurable;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.seqdoop.hadoop_bam.KeyIgnoringBAMOutputFormat;
import org.seqdoop.hadoop_bam.SAMRecordWritable;

import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.Base64;

/**
 * Takes a few Reads and will write them to a BAM file.
 * The Reads don't have to be sorted initially, the BAM file will be.
 * All the reads must fit into a single worker's memory, so this won't go well if you have too many.
 */
public class SmallBamWriter implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Takes a few Reads and will write them to a BAM file.
     * The Reads don't have to be sorted initially, the BAM file will be.
     * All the reads must fit into a single worker's memory, so this won't go well if you have too many.
     *
     * @param pipeline the pipeline to add this operation to.
     * @param reads  the reads to write (they don't need to be sorted).
     * @param header the header that corresponds to the reads.
     * @param destPath the GCS or local path to write to (must start with "gs://" if writing to GCS).
     * @param parquet whether to write out BAM or Parquet data (BDG AlignmentRecords); only applies when writing to Hadoop
     */
    public static void writeToFile(
            Pipeline pipeline, PCollection<GATKRead> reads, final SAMFileHeader header, final String destPath,
            final boolean parquet) {
        if ( BucketUtils.isHadoopUrl(destPath) ||
                pipeline.getRunner().getClass().equals(SparkPipelineRunner.class)) {
            writeToHadoop(pipeline, reads, header, destPath, parquet);
        } else {
            PCollectionView<Iterable<GATKRead>> iterableView =
                    reads.apply(View.<GATKRead>asIterable());

            PCollection<String> dummy = pipeline.apply("output file name", Create.<String>of(destPath));

            dummy.apply(ParDo.named("save to BAM file")
                            .withSideInputs(iterableView)
                            .of(new SaveToBAMFile(header, iterableView))
            );
        }
    }

    public static void writeToFile(Pipeline pipeline, PCollection<GATKRead> reads, final SAMFileHeader header, final String destPath) {
        writeToFile(pipeline, reads, header, destPath, false);
    }

    private static class SaveToBAMFile extends DoFn<String,Void> {
        private static final Logger logger = LogManager.getLogger(SaveToBAMFile.class);
        private static final long serialVersionUID = 1L;
        private final SAMFileHeader header;
        private final PCollectionView<Iterable<GATKRead>> iterableView;

        public SaveToBAMFile(SAMFileHeader header, PCollectionView<Iterable<GATKRead>> iterableView) {
            this.header = header;
            this.iterableView = iterableView;
        }

        @Override
        public void processElement(ProcessContext c) throws Exception {
            String dest = c.element();
            logger.info("Saving to " + dest);
            Iterable<GATKRead> reads = c.sideInput(iterableView);
            OutputStream outputStream = BucketUtils.createFile(dest, c.getPipelineOptions());
            try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, false, outputStream)) {
                for (GATKRead r : reads) {
                    final SAMRecord sr = r.convertToSAMRecord(header);
                    writer.addAlignment(sr);
                }
            }
        }
    }

    private static void writeToHadoop(
            Pipeline pipeline, PCollection<GATKRead> reads, final SAMFileHeader header, final String destPath,
            final boolean parquet) {
        if (destPath.equals("/dev/null")) {
            return;
        }

        String headerString = Base64.getEncoder().encodeToString(SerializableUtils.serializeToByteArray(header));

        @SuppressWarnings("unchecked")
        Class<? extends FileOutputFormat<NullWritable, SAMRecordWritable>> outputFormatClass =
                (Class<? extends FileOutputFormat<NullWritable, SAMRecordWritable>>) (Class<?>) TemplatedKeyIgnoringBAMOutputFormat.class;
        @SuppressWarnings("unchecked")
        HadoopIO.Write.Bound<NullWritable, SAMRecordWritable> write = HadoopIO.Write.to(destPath,
                outputFormatClass, NullWritable.class, SAMRecordWritable.class)
                .withConfigurationProperty(TemplatedKeyIgnoringBAMOutputFormat.SAM_HEADER_PROPERTY_NAME, headerString);

        PCollection<KV<NullWritable, SAMRecordWritable>> samReads =
                reads.apply(ParDo.of(new DoFn<GATKRead, KV<NullWritable, SAMRecordWritable>>() {
            private static final long serialVersionUID = 1L;

            @Override
            public void processElement(ProcessContext c) throws Exception {
                SAMRecord samRecord = c.element().convertToSAMRecord(header);
                SAMRecordWritable samRecordWritable = new SAMRecordWritable();
                samRecordWritable.set(samRecord);
                c.output(KV.of(NullWritable.get(), samRecordWritable));
            }
        })).setCoder(KvCoder.of(WritableCoder.of(NullWritable.class), WritableCoder.of(SAMRecordWritable.class)));

        // write as a single (unsharded) file
        samReads.apply(write.withoutSharding());
    }

    public static class TemplatedKeyIgnoringBAMOutputFormat<K> extends KeyIgnoringBAMOutputFormat<K>
        implements Configurable, ShardNameTemplateAware {

        public static final String SAM_HEADER_PROPERTY_NAME = "sam.header";

        private Configuration conf;

        @Override
        public void setConf(Configuration conf) {
            this.conf = conf;
            if (conf != null) {
                String headerString = conf.get(SAM_HEADER_PROPERTY_NAME);
                if (headerString == null) {
                    throw new IllegalStateException("SAM file header has not been set");
                }
                byte[] headerBytes = Base64.getDecoder().decode(headerString);
                setSAMHeader((SAMFileHeader) SerializableUtils.deserializeFromByteArray(headerBytes, "SAMFileHeader"));
            }
        }

        @Override
        public Configuration getConf() {
            return getConf();
        }

        @Override
        public void checkOutputSpecs(JobContext job) {
            // overwrite old files if present (consistent with dataflow behavior)
        }

        @Override
        public Path getDefaultWorkFile(TaskAttemptContext context,
                String extension) throws IOException {
            // Use ShardNameTemplateHelper to construct the output file name specified in HadoopIO.Write.to()
            // above, along with any sharding and suffix information (if specified).
            return ShardNameTemplateHelper.getDefaultWorkFile(this, context);
        }
    }

}
