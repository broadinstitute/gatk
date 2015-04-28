package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.tools.FlagStat;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class FlagStatDataflowUnitTest {

    @Test(groups = "dataflow")
    public void testFlagStatDataflowTransform(){
        File bam = new File(BaseTest.publicTestDir, "org/broadinstitute/hellbender/tools/dataflow/pipelines/FlagStatDataflow/flag_stat.bam");
        Pipeline p = TestPipeline.create();
        DataflowWorkarounds.registerGenomicsCoders(p);
        List<SimpleInterval> intervals = Lists.newArrayList(new SimpleInterval("chr1", 1, 101),
                new SimpleInterval("chr2", 1, 101),
                new SimpleInterval("chr3", 1, 101),
                new SimpleInterval("chr4", 1, 101),
                new SimpleInterval("chr5", 1, 101),
                new SimpleInterval("chr6", 1, 202),
                new SimpleInterval("chr7", 1, 202),
                new SimpleInterval("chr8", 1, 202));
        PCollection<Read> preads = DataflowUtils.getReadsFromLocalBams(p, intervals, Lists.newArrayList(bam));
        PCollection<FlagStat.FlagStatus> presult = preads.apply(new FlagStatusDataflowTransform());


        ReadsDataSource nonDataflow = new ReadsDataSource(bam);
        nonDataflow.setIntervalsForTraversal(intervals);
        FlagStat.FlagStatus expectedNonDataflow = new FlagStat.FlagStatus();
        for (SAMRecord sam: nonDataflow){
            expectedNonDataflow.add(sam);
        }

        DataflowAssert.thatSingleton(presult).isEqualTo(expectedNonDataflow);
        p.run();


    }
}
