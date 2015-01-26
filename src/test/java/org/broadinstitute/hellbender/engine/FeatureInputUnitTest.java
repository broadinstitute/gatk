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

public class FeatureInputUnitTest extends BaseTest {

    @DataProvider(name = "InvalidFeatureArgumentValuesDataProvider")
    public Object[][] getInvalidFeatureArgumentValues() {
        return new Object[][] {
                { "name:file:file2" },
                { "name:" },
                { ":file" },
                { ":" },
                { "::" },
                { "" }
        };
    }

    @Test(dataProvider = "InvalidFeatureArgumentValuesDataProvider", expectedExceptions = UserException.BadArgumentValue.class)
    public void testInvalidFeatureArgumentValue( final String invalidFeatureArgumentValue ) {
        FeatureInput<Feature> featureInput = new FeatureInput<>(invalidFeatureArgumentValue);
    }

    @Test
    public void testNoFeatureNameSpecified() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("myFile");

        Assert.assertEquals(featureInput.getFeatureFile(), new File("myFile"), "Wrong File in FeatureInput");
        // Name should default to the absolute path of the File when no name is specified
        Assert.assertEquals(featureInput.getName(), new File("myFile").getAbsolutePath(), "Wrong default name in FeatureInput");
    }

    @Test
    public void testFeatureNameSpecified() {
        FeatureInput<Feature> featureInput = new FeatureInput<>("myName:myFile");

        Assert.assertEquals(featureInput.getFeatureFile(), new File("myFile"), "Wrong File in FeatureInput");
        Assert.assertEquals(featureInput.getName(), "myName", "Wrong name in FeatureInput");
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
}
