package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.Level;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Map;

public class SparkCommandLineArgumentCollectionTest {

    private static final String prop1 = "spark.config";
    private static final String value1 = "something";
    private static final String prop2 = "spark.someother.value";
    private static final String value2 = "somethingelse";


    @Test(groups = "spark")
    public void testGetSparkProperties(){
        final SparkCommandLineArgumentCollection sparkArgumentCollection = new SparkCommandLineArgumentCollection();

        sparkArgumentCollection.sparkProperties.addAll(Arrays.asList(prop1+ '=' +value1, prop2 + '=' + value2));
        final Map<String, String> sparkProperties = sparkArgumentCollection.getSparkProperties();
        Assert.assertEquals(sparkProperties.size(),2);
        Assert.assertEquals(sparkProperties.get(prop1), value1);
        Assert.assertEquals(sparkProperties.get(prop2), value2);
    }

    @DataProvider(name="badSplits")
    public Object[][] badSplits(){
        return new Object[][] {
                {"spark"},
                {"="},
                {"spark="},
                {"=spark"},
                {"spark=spark=spark"}
        };
    }

    @Test(dataProvider = "badSplits", expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testBadProperties(String property){
        final SparkCommandLineArgumentCollection sparkArgumentCollection = new SparkCommandLineArgumentCollection();
        sparkArgumentCollection.sparkProperties.add(property);
        sparkArgumentCollection.getSparkProperties();
    }

    @DataProvider(name="toolLogLevels")
    public Object[][] toolLogLevels(){
        return new Object[][] {
                {Log.LogLevel.DEBUG, Level.DEBUG.name()},
                {Log.LogLevel.INFO, Level.INFO.name()},
                {Log.LogLevel.WARNING, Level.WARN.name()},
                {Log.LogLevel.ERROR, Level.ERROR.name()},
        };
    }

    @Test(dataProvider = "toolLogLevels")
    public void testSparkVerbosityWithValidToolLogLevels(Log.LogLevel toolVerbosity, String expectedSparkLevel){
        final SparkCommandLineArgumentCollection sparkArgumentCollection = new SparkCommandLineArgumentCollection();
        final String level = sparkArgumentCollection.getSparkVerbosity(toolVerbosity);
        Assert.assertEquals(level, expectedSparkLevel);
    }

    @DataProvider(name="sparkLogLevels")
    public Object[][] sparkLogLevels(){
        return new Object[][] {
                {Level.ALL.name()},
                {Level.DEBUG.name()},
                {Level.INFO.name()},
                {Level.WARN.name()},
                {Level.ERROR.name()},
                {Level.FATAL.name()},
                {Level.TRACE.name()},
                {Level.OFF.name()}
        };
    }

    @Test(dataProvider = "sparkLogLevels")
    public void test(String sparkLevel){
        final SparkCommandLineArgumentCollection sparkArgumentCollection = new SparkCommandLineArgumentCollection();
        sparkArgumentCollection.setSparkVerbosity(sparkLevel);
        final String level = sparkArgumentCollection.getSparkVerbosity(Log.LogLevel.INFO);
        Assert.assertEquals(level, sparkLevel);
    }

}
