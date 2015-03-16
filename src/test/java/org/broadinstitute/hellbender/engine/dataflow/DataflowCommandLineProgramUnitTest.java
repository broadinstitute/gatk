package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.testng.annotations.Test;


public class DataflowCommandLineProgramUnitTest {

    @Test(expectedExceptions = UserException.class)
    public void testUserExceptionUnwrapping(){
        Pipeline p = TestPipeline.create();
        PCollection<String> pstrings = p.apply(Create.of(Lists.newArrayList("Some", "Values")));
        pstrings.apply(DataflowUtils.throwException(new UserException("fail")));
        DataflowCommandLineProgram.runPipeline(p);
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testOtherExceptionsNotConvertedToUserExceptions(){
        Pipeline p = TestPipeline.create();
        PCollection<String> pstrings = p.apply(Create.of(Lists.newArrayList("Some","Values")));
        pstrings.apply(DataflowUtils.throwException(new GATKException("fail")));
        DataflowCommandLineProgram.runPipeline(p);

    }

}
