package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.cloudera.dataflow.hadoop.HadoopIO;
import com.cloudera.dataflow.hadoop.WritableCoder;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.base.Joiner;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.hadoop.io.LongWritable;
import org.broadinstitute.hellbender.engine.dataflow.coders.VariantCoder;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.seqdoop.hadoop_bam.VCFInputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.io.Serializable;
import java.util.List;

/**
 * Loads variants into PCollections using a Hadoop filesystem.
 */
public final class VariantsHadoopSource implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Reads variants from the specified files.
     * @param commaSeparatedPaths one or more Hadoop filesystem paths, separated by commas
     * @param pipeline the pipeline
     * @return a PCollection of variants found in the file.
     */
    public static PCollection<Variant> getAllVariants(String commaSeparatedPaths, Pipeline pipeline) {
        PCollection<KV<LongWritable, VariantContextWritable>> input = pipeline.apply(
                HadoopIO.Read.from(commaSeparatedPaths, VCFInputFormat.class, LongWritable.class, VariantContextWritable.class));
        input.setCoder(KvCoder.of(WritableCoder.of(LongWritable.class), WritableCoder.of(VariantContextWritable.class)));
        return input.apply(ParDo.of(new ConvertToVariant())).setCoder(new VariantCoder());
    }

    /**
     * Reads variants from the specified files.
     * @param paths one or more Hadoop filesystem paths
     * @param pipeline the pipeline
     * @return a PCollection of variants found in the file.
     */
    public static PCollection<Variant> getAllVariants(List<String> paths, Pipeline pipeline) {
        return getAllVariants(Joiner.on(",").join(paths), pipeline);
    }

    static class ConvertToVariant extends DoFn<KV<LongWritable, VariantContextWritable>, Variant> {
        private static final long serialVersionUID = 1L;
        @Override
        public void processElement(ProcessContext c) throws Exception {
            VariantContext variantContext = c.element().getValue().get();
            c.output(new VariantContextVariantAdapter(variantContext));
        }
    }
}
