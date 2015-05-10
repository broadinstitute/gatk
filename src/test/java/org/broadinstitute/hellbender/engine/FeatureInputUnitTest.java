package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.table.TableFeature;
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
                { "name,key=value,key=value:file" }      //duplicate key
        };
    }

    @Test(dataProvider = "InvalidFeatureArgumentValuesDataProvider", expectedExceptions = UserException.BadArgumentValue.class)
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

        Assert.assertEquals(featureInput.getFeatureFile(), new File(validFileOnlyFeatureArgumentValue), "Wrong File in FeatureInput");
        // Name should default to the absolute path of the File when no name is specified
        Assert.assertEquals(featureInput.getName(), new File(validFileOnlyFeatureArgumentValue).getAbsolutePath(), "Wrong default name in FeatureInput");
    }

    @Test
    public void testFeatureNameSpecified() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("myName:myFile");

        Assert.assertEquals(featureInput.getFeatureFile(), new File("myFile"), "Wrong File in FeatureInput");
        Assert.assertEquals(featureInput.getName(), "myName", "Wrong name in FeatureInput");
    }

    @Test
    public void testNullOKAsFeatureName() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("null:myFile");

        Assert.assertEquals(featureInput.getFeatureFile(), new File("myFile"), "Wrong File in FeatureInput");
        Assert.assertEquals(featureInput.getName(), "null", "Wrong name in FeatureInput");
    }

    @Test
    public void testNullOKAsFileName() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("myName:null");

        Assert.assertEquals(featureInput.getFeatureFile(), new File("null"), "Wrong File in FeatureInput");
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
        Assert.assertEquals(featureInput.getFeatureFile(), new File("myFile"));
    }

    @Test
    public void testFeatureKeyValuePairSpecified() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("myName,key1=value1:myFile");

        Assert.assertEquals(featureInput.getAttribute("key1"), "value1", "wrong attribute value for key1");
        Assert.assertEquals(featureInput.getAttribute("key2"), null, "wrong attribute value for key2 (not present)");

        Assert.assertEquals(featureInput.getName(), "myName");
        Assert.assertEquals(featureInput.getFeatureFile(), new File("myFile"));
    }

    @Test
    public void testFeatureKeyValuePairsSpecifiedSameValue() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("myName,key1=value,key2=value:myFile");

        Assert.assertEquals(featureInput.getAttribute("key1"), "value", "wrong attribute value for key1");
        Assert.assertEquals(featureInput.getAttribute("key2"), "value", "wrong attribute value for key2");

        Assert.assertEquals(featureInput.getName(), "myName");
        Assert.assertEquals(featureInput.getFeatureFile(), new File("myFile"));
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
    public void testGetFeatureType() {
        FeatureInput<Feature> featureFeatureInput = new FeatureInput<>("file1");
        FeatureInput<VariantContext> variantContextFeatureInput = new FeatureInput<>("file2");
        FeatureInput<TableFeature> tableFeatureInput = new FeatureInput<>("file3");

        featureFeatureInput.setFeatureType(Feature.class);
        variantContextFeatureInput.setFeatureType(VariantContext.class);
        tableFeatureInput.setFeatureType(TableFeature.class);

        Assert.assertEquals(featureFeatureInput.getFeatureType(), Feature.class, "Wrong Feature type for FeatureInput<Feature>");
        Assert.assertEquals(variantContextFeatureInput.getFeatureType(), VariantContext.class, "Wrong Feature type for FeatureInput<VariantContext>");
        Assert.assertEquals(tableFeatureInput.getFeatureType(), TableFeature.class, "Wrong Feature type for FeatureInput<TableFeature>");
    }

    @Test
    public void testToString() {
        final FeatureInput<Feature> namelessFeatureInput = new FeatureInput<>("file1");
        final FeatureInput<Feature> namedFeatureInput = new FeatureInput<>("name:file1");

        Assert.assertEquals(namelessFeatureInput.toString(), new File("file1").getAbsolutePath(), "String representation of nameless FeatureInput incorrect");
        Assert.assertEquals(namedFeatureInput.toString(), "name:" + new File("file1").getAbsolutePath(), "String representation of named FeatureInput incorrect");
    }
}
