package org.broadinstitute.hellbender.engine.spark;

import org.apache.hadoop.test.MockitoMaker;
import org.apache.spark.SparkContext;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.tools.ant.taskdefs.Java;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.TestSparkProgramGroup;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;


public class SparkCommandLineProgramUnitTest extends BaseTest {

    @CommandLineProgramProperties(summary = "test", programGroup = TestSparkProgramGroup.class, oneLineSummary = "line")
    private static class TestSparkCommandLineProgram extends SparkCommandLineProgram {
        private static final long serialVersionUID = 1L;
        @Override
        protected void runPipeline(JavaSparkContext ctx) {
        }
    };

    // Because the exact values here are system dependant in local mode, this tests are just sanity tests
    @Test (dataProvider = "testMemoryPerCoreData")
    public void testMemoryPerCore(JavaSparkContext context, double expected) {
        TestSparkCommandLineProgram prog = new TestSparkCommandLineProgram();
        prog.ctx = context;
        Assert.assertNotEquals(prog.getMemoryPerCore(),0);
        Assert.assertEquals(prog.getMemoryPerCore()>0,true);
    }



    @Test (dataProvider = "testMemoryPerCoreData")
    public void testCoresPerExecutor(JavaSparkContext context, double expected) {
        TestSparkCommandLineProgram prog = new TestSparkCommandLineProgram();
        prog.ctx = context;
        Assert.assertNotNull(prog.getExecutorCores());
        Assert.assertNotEquals(prog.getExecutorCores(), 0);
        Assert.assertEquals(prog.getExecutorCores()>0, true);
    }

    @DataProvider (name = "testMemoryPerCoreData")
    public Object[][] testMemoryPerCoreData() {
        return new Object[][] {
                {makeSparkContext("1", "1G", "3", "4", "12G", true), 3072.0}
        };
    }

    private JavaSparkContext makeSparkContext(String driverCores, String driverMemory,
                                  String numExecutors, String executorCores, String executorMemory, boolean local) {
        Map<String, String> properties = new HashMap<>();
        properties.put("spark.driver.memory", driverMemory);
        properties.put("spark.driver.cores", driverCores);
        properties.put("spark.executor.instances", numExecutors);
        properties.put("spark.executor.memory", executorMemory);
        properties.put("spark.executor.cores", executorCores);
        if (!local) {
            properties.put("spark.master", "spark://127.0.0.1:7077");
        }
        return SparkContextFactory.getTestSparkContext(properties);

    }
}
