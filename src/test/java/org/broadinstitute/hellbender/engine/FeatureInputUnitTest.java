package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class FeatureInputUnitTest extends BaseTest {

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
                { "name,key=value,key=value:file" },      //duplicate key
                { "name:name:gendb://mydb" }
        };
    }

    @Test(dataProvider = "InvalidFeatureArgumentValuesDataProvider", expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testInvalidFeatureArgumentValue( final String invalidFeatureArgumentValue ) {
        FeatureInput<Feature> featureInput = new FeatureInput<>(invalidFeatureArgumentValue);
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
        FeatureInput<Feature> featureInput = new FeatureInput<>(validFileOnlyFeatureArgumentValue);   //"myName,key1=value,myFile"

        Assert.assertEquals(featureInput.getFeaturePath(), validFileOnlyFeatureArgumentValue, "Wrong File in FeatureInput");
        // Name should default to the absolute path of the File when no name is specified
        Assert.assertEquals(featureInput.getName(), new File(validFileOnlyFeatureArgumentValue).getAbsolutePath(), "Wrong default name in FeatureInput");
    }

    @DataProvider(name = "GenDbPathAndNameData")
    public Object[][] genDbPathAndNameData() {
        return new Object[][] {
                // input String, expected Feature path, expected logical name
                {"gendb://myJsons", "gendb://myJsons", "gendb://" + new File("myJsons").getAbsolutePath()},
                {"myname:gendb://myJsons", "gendb://myJsons", "myname"},
                {"myname,key1=value1:gendb://myJsons", "gendb://myJsons", "myname"},
                {"myname//:gendb://myJsons", "gendb://myJsons", "myname//"},
                {"myname:gendb://", "gendb://", "myname"}
        };
    }

    @Test(dataProvider = "GenDbPathAndNameData")
    public void testGenDbPathAndName( final String inputString, final String expectedFeaturePath, final String expectedLogicalName ) {
        final FeatureInput<VariantContext> gendbInput = new FeatureInput<>(inputString);

        Assert.assertEquals(gendbInput.getFeaturePath(), expectedFeaturePath, "wrong featurePath");
        Assert.assertEquals(gendbInput.getName(), expectedLogicalName, "wrong logical name");
    }

    @Test
    public void testFeatureNameSpecified() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("myName:myFile");

        Assert.assertEquals(featureInput.getFeaturePath(), "myFile", "Wrong File in FeatureInput");
        Assert.assertEquals(featureInput.getName(), "myName", "Wrong name in FeatureInput");
    }

    @Test
    public void testNullOKAsFeatureName() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("null:myFile");

        Assert.assertEquals(featureInput.getFeaturePath(), "myFile", "Wrong File in FeatureInput");
        Assert.assertEquals(featureInput.getName(), "null", "Wrong name in FeatureInput");
    }

    @Test
    public void testNullOKAsFileName() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("myName:null");

        Assert.assertEquals(featureInput.getFeaturePath(), "null", "Wrong File in FeatureInput");
        Assert.assertEquals(featureInput.getName(), "myName", "Wrong name in FeatureInput");
    }


    @Test
    public void testFeatureKeyValuePairsSpecified() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("myName,key1=value1,key2=value2,null=null:myFile");

        Assert.assertEquals(featureInput.getAttribute("key1"), "value1", "wrong attribute value for key1");
        Assert.assertEquals(featureInput.getAttribute("key2"), "value2", "wrong attribute value for key2");
        Assert.assertEquals(featureInput.getAttribute("null"), "null", "wrong attribute value for key \"null\"");
        Assert.assertEquals(featureInput.getAttribute("key3"), null, "wrong attribute value for key3 (not present)");

        Assert.assertEquals(featureInput.getName(), "myName");
        Assert.assertEquals(featureInput.getFeaturePath(), "myFile");
    }

    @Test
    public void testFeatureKeyValuePairSpecified() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("myName,key1=value1:myFile");

        Assert.assertEquals(featureInput.getAttribute("key1"), "value1", "wrong attribute value for key1");
        Assert.assertEquals(featureInput.getAttribute("key2"), null, "wrong attribute value for key2 (not present)");

        Assert.assertEquals(featureInput.getName(), "myName");
        Assert.assertEquals(featureInput.getFeaturePath(), "myFile");
    }

    @Test
    public void testFeatureKeyValuePairsSpecifiedSameValue() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("myName,key1=value,key2=value:myFile");

        Assert.assertEquals(featureInput.getAttribute("key1"), "value", "wrong attribute value for key1");
        Assert.assertEquals(featureInput.getAttribute("key2"), "value", "wrong attribute value for key2");

        Assert.assertEquals(featureInput.getName(), "myName");
        Assert.assertEquals(featureInput.getFeaturePath(), "myFile");
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
    public void testFeatureValuesForNullKey(final String featureInputArgument ) {
        FeatureInput<Feature> featureInput = new FeatureInput<>(featureInputArgument);
        featureInput.getAttribute(null);
    }

    @Test
    public void testToString() {
        final FeatureInput<Feature> namelessFeatureInput = new FeatureInput<>("file1");
        final FeatureInput<Feature> namedFeatureInput = new FeatureInput<>("name:file1");
        final FeatureInput<Feature> namelessGenomicsDB = new FeatureInput<>("gendb://file1");
        final FeatureInput<Feature> namedGenomicsDB = new FeatureInput<>("name:gendb://file1");

        Assert.assertEquals(namelessFeatureInput.toString(), new File("file1").getAbsolutePath(), "String representation of nameless FeatureInput incorrect");
        Assert.assertEquals(namedFeatureInput.toString(), "name:" + new File("file1").getAbsolutePath(), "String representation of named FeatureInput incorrect");
        Assert.assertEquals(namelessGenomicsDB.toString(), "gendb://" + new File("file1").getAbsolutePath(), "String representation of nameless FeatureInput with genomicsDB path incorrect");
        Assert.assertEquals(namedGenomicsDB.toString(), "name:gendb://" + new File("file1").getAbsolutePath(), "String representation of named FeatureInput with genomicsDB path incorrect");
    }
}
