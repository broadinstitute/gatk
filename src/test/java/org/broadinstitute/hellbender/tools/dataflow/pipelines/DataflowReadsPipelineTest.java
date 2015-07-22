package org.broadinstitute.hellbender.tools.dataflow.pipelines;


import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.dataflow.transforms.CountBasesDataflowTransform;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;
import java.util.List;


public final class DataflowReadsPipelineTest {

    @CommandLineProgramProperties(summary ="runtime configurable pipeline for testing", oneLineSummary = "for testing pipeline functionality")
    private static class RuntimeConfigurablePipeline extends DataflowReadsPipeline {
        private static final long serialVersionUID = 1l;
        private final PTransformSAM<?> tool;
        private final List<ReadFilter> filters;
        private final List<ReadTransformer> transformers;

        public RuntimeConfigurablePipeline(PTransformSAM<?> tool, List<ReadFilter> filters, List<ReadTransformer> transformers ){
            this.tool = tool;
            this.filters = filters;
            this.transformers = transformers;
        }

        @Override
        protected ImmutableList<ReadFilter> getReadFilters( SAMFileHeader header ) {
            return ImmutableList.copyOf(filters);
        }

        @Override
        protected ImmutableList<ReadTransformer> getReadTransformers(){
            return ImmutableList.copyOf(transformers);
        }

        @Override
        protected PTransformSAM<?> getTool() {
            return tool;
        }
    }

    @SuppressWarnings("unchecked")
    @DataProvider(name = "configurations")
    public Object[][] configurations(){
        ReadFilter all = r -> true;
        ReadFilter none = r -> false;
        ReadTransformer replaceBasesWithA = r -> {
                r.setBases(new byte[]{(byte) 'A'});
                return r;
        };

        return new Object[][] {
                new Object[] { Collections.<ReadFilter>emptyList(), Collections.<ReadTransformer>emptyList(), 707},
                new Object[]{ Lists.newArrayList(all), Collections.emptyList(), 707},
                new Object[] { Lists.newArrayList(none), Collections.emptyList(), 0},
                new Object[] { Collections.<ReadFilter>emptyList(), Lists.newArrayList(replaceBasesWithA), 7},
                new Object[] { Lists.newArrayList(all), Lists.newArrayList(replaceBasesWithA), 7},
                new Object[] { Lists.newArrayList(all), Lists.newArrayList(replaceBasesWithA.andThen(ReadTransformer.identity())), 7},
                new Object[] { Lists.newArrayList(none.negate().and(all)), Lists.newArrayList(replaceBasesWithA.compose(ReadTransformer.identity())), 7},
                new Object[] { Lists.newArrayList(none.negate(),all, all), Lists.newArrayList(replaceBasesWithA, replaceBasesWithA, ReadTransformer.identity()), 7}
        };
    }

    @SuppressWarnings("unchecked")
    @Test(dataProvider = "configurations", groups = {"dataflow"})
    public void testTransformersAndFilters(List<ReadFilter> filters, List<ReadTransformer> transformers, long expectedCounts){
        RuntimeConfigurablePipeline rcp = new RuntimeConfigurablePipeline(new CountBasesDataflowTransform(), filters, transformers);
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);
        String bam = "src/test/resources/org/broadinstitute/hellbender/tools/count_reads_sorted.bam";
        PCollection<GATKRead> preads = DataflowUtils.getReadsFromLocalBams(p, Lists.newArrayList(new SimpleInterval("chr7", 1, 404)), Lists.newArrayList(new File(bam)));

        PCollection<?> presult = rcp.applyTransformsToPipeline(SamReaderFactory.makeDefault().getFileHeader(new File(bam)), preads);

        DataflowAssert.thatSingleton((PCollection<Long>) presult).isEqualTo(expectedCounts);

        p.run();
    }
}
