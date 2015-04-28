package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.appengine.repackaged.com.google.common.collect.Lists;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;
import java.util.stream.Collectors;

public class CountBasesDataflowUnitTest  {

    @DataProvider(name = "sams")
    public Object[][] sams(){

        return new Object[][]{
                {Lists.newArrayList(ArtificialSAMUtils.createRandomRead(100)), 100,},
                {Lists.newArrayList(ArtificialSAMUtils.createRandomRead(50), ArtificialSAMUtils.createRandomRead(100), ArtificialSAMUtils.createRandomRead(1000)), 1150}
        };
    }



    @Test(dataProvider = "sams", groups= {"dataflow"})
    public void countBasesTest(List<SAMRecord> sams, long expected){
        List<Read> reads = sams.stream().map(ReadConverter::makeRead).collect(Collectors.toList());
        Pipeline p = TestPipeline.create();
        DataflowWorkarounds.registerGenomicsCoders(p);
        PCollection<Read> preads = p.apply(Create.of(reads));
        PTransformSAM<Long> transform = new CountBasesDataflowTransform();
        transform.setHeaderString(ArtificialSAMUtils.createArtificialSamHeader().toString());
        PCollection<Long> presult = preads.apply(transform);

        DataflowAssert.thatSingleton(presult).isEqualTo(expected);

        p.run();

    }

}