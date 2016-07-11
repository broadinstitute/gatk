package org.broadinstitute.hellbender.engine;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class TaggedInputFileArgumentUnitTest extends BaseTest {

    @CommandLineProgramProperties(
            summary = "test",
            oneLineSummary = "",
            programGroup = TestProgramGroup.class
    )
    private static final class TestReadWalker extends ReadWalker{
        TaggedInputFileArgument tumor;
        TaggedInputFileArgument normal;
        TaggedInputFileArgument met;

        @Override
        public void onTraversalStart() {
            tumor  = readArguments.getInputsBySymbolicName().get("tumor");
            normal = readArguments.getInputsBySymbolicName().get("normal");
            met = readArguments.getInputsBySymbolicName().get("met");
        }

        @Override
        public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
             //do nothing
        }
    }

    @Test
    public void testSymbolicNames() throws Exception {
        final TestReadWalker tool = new TestReadWalker();
        final CommandLineParser clp = new CommandLineParser(tool);
        final File bamFileT = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
        final File bamFileN = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test2.bam");
        final String[] args = {
                "-I", "tumor:" + bamFileT.getCanonicalPath(),
                "-I", "normal:" + bamFileN.getCanonicalPath(),
        };
        clp.parseArguments(System.out, args);
        tool.onStartup();
        Assert.assertNull(tool.tumor, "tumor");
        Assert.assertNull(tool.normal, "normal");
        Assert.assertNull(tool.met, "met");
        tool.doWork();
        Assert.assertNotNull(tool.tumor, "tumor");
        Assert.assertEquals(tool.tumor.getFile(), bamFileT);
        Assert.assertNotNull(tool.normal, "normal");
        Assert.assertEquals(tool.normal.getFile(), bamFileN);
        Assert.assertNull(tool.met, "met");   //this one was not provided
        tool.onShutdown();
    }

    @CommandLineProgramProperties(summary = "",
            oneLineSummary = "",
            programGroup = TestProgramGroup.class)
    private static final class TestGATKSparkTool extends GATKSparkTool{
        private static final long serialVersionUID = -7467286444549225596L;
        TaggedInputFileArgument tumor;
        TaggedInputFileArgument normal;

        @Override
        protected void runTool( JavaSparkContext ctx ){
            tumor = readArguments.getInputsBySymbolicName().get("tumor");
            normal = readArguments.getInputsBySymbolicName().get("normal");
        }
    }

    @Test
    public void testSymbolicNamesInSpark() throws Exception {
        final TestGATKSparkTool tool = new TestGATKSparkTool();
        final CommandLineParser clp = new CommandLineParser(tool);
        final File bamFileT = new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");
        final String[] args = {
                "-I", "tumor:" + bamFileT.getCanonicalPath(),
        };
        //Note: only 1 input file is supported for Spark tools
        clp.parseArguments(System.out, args);
        Assert.assertNull(tool.tumor, "tumor");
        Assert.assertNull(tool.normal, "normal");   //this one was not provided
        tool.runTool();
        Assert.assertNotNull(tool.tumor, "tumor");
        Assert.assertEquals(tool.tumor.getFile(), bamFileT);
        Assert.assertNull(tool.normal, "normal");   //this one was not provided
    }

    @DataProvider(name = "InvalidFeatureArgumentValuesDataProvider")
    public Object[][] getInvalidFeatureArgumentValues() {
        return new Object[][] {
                { "name:file:file2" },
                { "name:" },
                { ":file" },
                { ",:file" },
                { "name,key=value=fred:file" },
                { "name,:file"},
                { ",key=value:file"},
                {"name,key:file"},
                {"name,key=:file"},
                {"name,=value:file"},
                {"name,=:file"},
                { ":" },
                { ",:" },
                { "::" },
                { "" },
                { "name,key=value1,key=value2:file" },   //duplicate key
                { "name,key=value,key=value:file" }      //duplicate key
        };
    }

    @Test(dataProvider = "InvalidFeatureArgumentValuesDataProvider", expectedExceptions = UserException.BadArgumentValue.class)
    public void testInvalidFeatureArgumentValue( final String invalidFeatureArgumentValue ) {
        TaggedInputFileArgument taggedInputFileArgument = new TaggedInputFileArgument(invalidFeatureArgumentValue);
    }

    @DataProvider(name = "ValidFileOnlyFeatureArgumentValuesDataProvider")
    public Object[][] getValidFileOnlyFeatureArgumentValues() {
        return new Object[][] {
                {"myFile"},
                {"myName,key1=value,myFile"},     //allowed - all of this is treated as a file name
                {"=myFile"},                      //allowed - all of this is treated as a file name
                {",myFile"},                    //allowed - all of this is treated as a file name
                {"=,myFile"},                    //allowed - all of this is treated as a file name
                {"key1=value,myFile"}             //allowed - all of this is treated as a file name
        };
    }
    @Test(dataProvider = "ValidFileOnlyFeatureArgumentValuesDataProvider")
    public void testNoFeatureNameSpecified(final String validFileOnlyFeatureArgumentValue) {
        TaggedInputFileArgument taggedInputFileArgument = new TaggedInputFileArgument(validFileOnlyFeatureArgumentValue);   //"myName,key1=value,myFile"

        Assert.assertEquals(taggedInputFileArgument.getFile(), new File(validFileOnlyFeatureArgumentValue), "Wrong File in TaggedInputFileArgument");
        // Name should default to the absolute path of the File when no name is specified
        Assert.assertEquals(taggedInputFileArgument.getName(), new File(validFileOnlyFeatureArgumentValue).getAbsolutePath(), "Wrong default name in TaggedInputFileArgument");
    }

    @DataProvider(name = "ValidPairsForComparisonDataProvider")
    public Object[][] getValidPairsForComparison() {
        final String n1 = "myFile";
        final String n2 = "myName,key1=value,myFile";
        final String n3 = "myName,key1=value1,key2=value2:myFile";
        final String n3_noKV = "myName:myFile";
        final String n3_moreKV = "myName,key1=value1,key2=value2,key3=value3:myFile";
        return new Object[][] {
                {n1, n2, false},
                {n2, n3, false},
                {n1, n3, false},
                {n3, n3_noKV, true},
                {n3, n3_moreKV, true},
                {n3_noKV, n3_noKV, true},
        };
    }
    @Test(dataProvider = "ValidPairsForComparisonDataProvider")
    public void testEquality(final String string1, final String string2, final boolean expectedEquals) {
        final TaggedInputFileArgument tifa1 = new TaggedInputFileArgument(string1);
        final TaggedInputFileArgument tifa2 = new TaggedInputFileArgument(string2);
        Assert.assertNotEquals(tifa1, null);
        Assert.assertNotEquals(null, tifa1);
        Assert.assertNotEquals(tifa2, null);
        Assert.assertNotEquals(null, tifa2);
        Assert.assertEquals(tifa1, tifa1);
        Assert.assertEquals(tifa2, tifa2);
        Assert.assertNotEquals(tifa1, "fred");
        Assert.assertNotEquals(tifa2, "fred");
        Assert.assertNotEquals("fred", tifa1);
        Assert.assertNotEquals("fred", tifa2);

        Assert.assertEquals(tifa1.hashCode(), tifa1.hashCode());
        Assert.assertEquals(tifa2.hashCode(), tifa2.hashCode());

        Assert.assertEquals(tifa1.equals(tifa2), expectedEquals);
        Assert.assertEquals(tifa2.equals(tifa1), expectedEquals);

        final TaggedInputFileArgument tifa1Copy = new TaggedInputFileArgument(tifa1.getName(), tifa1.getAttributes(), tifa1.getFile());
        final TaggedInputFileArgument tifa2Copy = new TaggedInputFileArgument(tifa2.getName(), tifa2.getAttributes(), tifa2.getFile());
        Assert.assertEquals(tifa1, tifa1Copy);
        Assert.assertEquals(tifa1Copy, tifa1);
        Assert.assertEquals(tifa2, tifa2Copy);
        Assert.assertEquals(tifa2Copy, tifa2);
    }

    @Test
    public void testFeatureNameSpecified() {
        TaggedInputFileArgument taggedInputFileArgument = new TaggedInputFileArgument("myName:myFile");

        Assert.assertEquals(taggedInputFileArgument.getFile(), new File("myFile"), "Wrong File in TaggedInputFileArgument");
        Assert.assertEquals(taggedInputFileArgument.getName(), "myName", "Wrong name in TaggedInputFileArgument");
    }

    @Test
    public void testNullOKAsFeatureName() {
        TaggedInputFileArgument taggedInputFileArgument = new TaggedInputFileArgument("null:myFile");

        Assert.assertEquals(taggedInputFileArgument.getFile(), new File("myFile"), "Wrong File in TaggedInputFileArgument");
        Assert.assertEquals(taggedInputFileArgument.getName(), "null", "Wrong name in TaggedInputFileArgument");
    }

    @Test
    public void testNullOKAsFileName() {
        TaggedInputFileArgument taggedInputFileArgument = new TaggedInputFileArgument("myName:null");

        Assert.assertEquals(taggedInputFileArgument.getFile(), new File("null"), "Wrong File in TaggedInputFileArgument");
        Assert.assertEquals(taggedInputFileArgument.getName(), "myName", "Wrong name in TaggedInputFileArgument");
    }


    @Test
    public void testFeatureKeyValuePairsSpecified() {
        TaggedInputFileArgument taggedInputFileArgument = new TaggedInputFileArgument("myName,key1=value1,key2=value2,null=null:myFile");

        Assert.assertEquals(taggedInputFileArgument.getAttribute("key1"), "value1", "wrong attribute value for key1");
        Assert.assertEquals(taggedInputFileArgument.getAttribute("key2"), "value2", "wrong attribute value for key2");
        Assert.assertEquals(taggedInputFileArgument.getAttribute("null"), "null", "wrong attribute value for key \"null\"");
        Assert.assertEquals(taggedInputFileArgument.getAttribute("key3"), null, "wrong attribute value for key3 (not present)");

        Assert.assertEquals(taggedInputFileArgument.getName(), "myName");
        Assert.assertEquals(taggedInputFileArgument.getFile(), new File("myFile"));
    }

    @Test
    public void testFeatureKeyValuePairSpecified() {
        TaggedInputFileArgument taggedInputFileArgument = new TaggedInputFileArgument("myName,key1=value1:myFile");

        Assert.assertEquals(taggedInputFileArgument.getAttribute("key1"), "value1", "wrong attribute value for key1");
        Assert.assertEquals(taggedInputFileArgument.getAttribute("key2"), null, "wrong attribute value for key2 (not present)");

        Assert.assertEquals(taggedInputFileArgument.getName(), "myName");
        Assert.assertEquals(taggedInputFileArgument.getFile(), new File("myFile"));
    }

    @Test
    public void testFeatureKeyValuePairsSpecifiedSameValue() {
        TaggedInputFileArgument taggedInputFileArgument = new TaggedInputFileArgument("myName,key1=value,key2=value:myFile");

        Assert.assertEquals(taggedInputFileArgument.getAttribute("key1"), "value", "wrong attribute value for key1");
        Assert.assertEquals(taggedInputFileArgument.getAttribute("key2"), "value", "wrong attribute value for key2");

        Assert.assertEquals(taggedInputFileArgument.getName(), "myName");
        Assert.assertEquals(taggedInputFileArgument.getFile(), new File("myFile"));
    }

    @DataProvider(name = "KeyValuesDataProviderForTestingNull")
    public Object[][] getKeyValuesDataProviderForTestingNull() {
        return new Object[][] {
                { "myName,key1=value1,key2=value2:myFile" },
                { "myName,null=value:myFile"},
                { "myName:myFile" },
                { "myFile" },
                { "null" },
                { "null:myFile" },
                { "null:null" },
        };
    }

    @Test(dataProvider = "KeyValuesDataProviderForTestingNull", expectedExceptions = IllegalArgumentException.class)
    public void testFeatureValuesForNullKey(final String taggedInputFileArgumentArgument ) {
        TaggedInputFileArgument taggedInputFileArgument = new TaggedInputFileArgument(taggedInputFileArgumentArgument);
        taggedInputFileArgument.getAttribute(null);
    }

    @Test
    public void testToString() {
        final TaggedInputFileArgument namelessTaggedInputFileArgument = new TaggedInputFileArgument("file1");
        final TaggedInputFileArgument namedTaggedInputFileArgument = new TaggedInputFileArgument("name:file1");

        Assert.assertEquals(namelessTaggedInputFileArgument.toString(), new File("file1").getAbsolutePath(), "String representation of nameless TaggedInputFileArgument incorrect");
        Assert.assertEquals(namedTaggedInputFileArgument.toString(), "name:" + new File("file1").getAbsolutePath(), "String representation of named TaggedInputFileArgument incorrect");
    }
}

