package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.appengine.repackaged.com.google.common.collect.Lists;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class PrintReadsDataflowUnitTest {

    @DataProvider(name = "sams")
    public Object[][] sams(){

        return new Object[][]{
                {Lists.newArrayList(ArtificialSAMUtils.createRandomRead(30))},
                {Lists.newArrayList(ArtificialSAMUtils.createRandomRead(50), ArtificialSAMUtils.createRandomRead(100), ArtificialSAMUtils.createRandomRead(1000))}
        };
    }



    @Test(dataProvider = "sams", groups= {"dataflow"})
    public void printReadsTransformTest(List<SAMRecord> sams){
        List<Read> reads = sams.stream().map(ReadConverter::makeRead).collect(Collectors.toList());
        Set<String> expected = sams.stream().
                map(r -> {r.setMateUnmappedFlag(true); return r;}) // this is a weird hack, mateUnmapped is undefined if there is no mate
                // but our standard if false, google genomics defaults to true instead, but both are ok...
        .map(SAMRecord::getSAMString).collect(Collectors.toSet());

        Pipeline p = TestPipeline.create();
        DataflowWorkarounds.registerGenomicsCoders(p);
        PCollection<Read> preads = p.apply(Create.of(reads));
        PTransformSAM<String> transform = new PrintReadsDataflowTransform();
        transform.setHeaderString(ArtificialSAMUtils.createArtificialSamHeader().toString());
        PCollection<String> presult = preads.apply(transform);

        DataflowAssert.that(presult);

        p.run();

    }

}