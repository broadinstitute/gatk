package org.broadinstitute.hellbender.utils.dataflow;


import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.OutputStream;
import java.io.Serializable;

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
     */
    public static void writeToFile(Pipeline pipeline, PCollection<GATKRead> reads, final SAMFileHeader header, final String destPath) {
        PCollectionView<Iterable<GATKRead>> iterableView =
                reads.apply(View.<GATKRead>asIterable());

        PCollection<String> dummy = pipeline.apply("output file name", Create.<String>of(destPath));

        dummy.apply(ParDo.named("save to BAM file")
                        .withSideInputs(iterableView)
                        .of(new SaveToBAMFile(header, iterableView))
        );

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
            try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, outputStream)) {
                for (GATKRead r : reads) {
                    final SAMRecord sr = r.convertToSAMRecord(header);
                    writer.addAlignment(sr);
                }
            }
        }
    }
}
