package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;

public class CountReadsDataflowUnitTest extends BaseTest
{

    @DataProvider(name="numberOfReads")
    public Object[][] numberOfReads(){
        return new Object[][]{
                {1},
                {2},
                {100},
                {1000}
        };
    }

    @Test(dataProvider = "numberOfReads", groups = "dataflow")
    public void testCountReadsTransform(long numberOfReads){
        List<Read> reads = Collections.nCopies((int)numberOfReads, new Read() );

        Pipeline p = TestPipeline.create();
        DataflowWorkarounds.registerGenomicsCoders(p);
        PTransform<PCollection<Read>, PCollection<Long>> countReads = new CountReadsDataflowTransform();
        PCollection<Long> presult = p.apply(Create.of(reads)).apply(countReads);

        DataflowAssert.thatSingleton(presult).isEqualTo(numberOfReads);
        p.run();
    }

}