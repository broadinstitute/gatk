package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;

public final class CountReadsTransformUnitTest extends BaseTest{

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
        List<GATKRead> reads = Collections.nCopies((int)numberOfReads, new GoogleGenomicsReadToGATKReadAdapter(new Read()) );

        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);
        PTransform<PCollection<GATKRead>, PCollection<Long>> countReads = new CountReadsDataflowTransform();
        PCollection<Long> presult = p.apply(Create.of(reads)).apply(countReads);

        DataflowAssert.thatSingleton(presult).isEqualTo(numberOfReads);
        p.run();
    }

}