package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;

public final class FeatureInputUnitTest extends GATKBaseTest {
    private static final String FEATURE_INPUT_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

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
    
    @DataProvider(name = "GcsPathAndNameData")
    public Object[][] gcsPathAndNameData() {
        return new Object[][] {
                // input String, expected Feature path, expected logical name
                {"gs://bucket/user/my.vcf", "gs://bucket/user/my.vcf", "gs://bucket/user/my.vcf"},
                {"myname:gs://bucket/user/my.vcf", "gs://bucket/user/my.vcf", "myname"},
                {"myname,key1=value1:gs://bucket/user/my.vcf", "gs://bucket/user/my.vcf", "myname"},
                {"myname//:gs://bucket/user/my.vcf", "gs://bucket/user/my.vcf", "myname//"}
        };
    }
    
    @DataProvider(name = "HdfsPathAndNameData")
    public Object[][] hdfsPathAndNameData() {
        return new Object[][] {
                // input String, expected Feature path, expected logical name
                {"hdfs://localhost/user/my.vcf", "hdfs://localhost/user/my.vcf", "hdfs://localhost:8020/user/my.vcf"},
                {"myname:hdfs://localhost/user/my.vcf", "hdfs://localhost/user/my.vcf", "myname"},
                {"myname,key1=value1:hdfs://localhost/user/my.vcf", "hdfs://localhost/user/my.vcf", "myname"},
                {"myname//:hdfs://localhost/user/my.vcf", "hdfs://localhost/user/my.vcf", "myname//"}
        };
    }

    @Test(dataProvider = "GenDbPathAndNameData")
    public void testGenDbPathAndName( final String inputString, final String expectedFeaturePath, final String expectedLogicalName ) {
        final FeatureInput<VariantContext> gendbInput = new FeatureInput<>(inputString);

        Assert.assertEquals(gendbInput.getFeaturePath(), expectedFeaturePath, "wrong featurePath");
        Assert.assertEquals(gendbInput.getName(), expectedLogicalName, "wrong logical name");
    }

    @Test(dataProvider = "GcsPathAndNameData", groups={"bucket"})
    public void testGcsPathAndName( final String inputString, final String expectedFeaturePath, final String expectedLogicalName ) {
        final FeatureInput<VariantContext> gcsInput = new FeatureInput<>(inputString);

        Assert.assertEquals(gcsInput.getFeaturePath(), expectedFeaturePath, "wrong featurePath");
        Assert.assertEquals(gcsInput.getName(), expectedLogicalName, "wrong logical name");
    }

    @Test(dataProvider = "HdfsPathAndNameData")
    public void testHdfsPathAndName( final String inputString, final String expectedFeaturePath, final String expectedLogicalName ) {
        final FeatureInput<VariantContext> hdfsInput = new FeatureInput<>(inputString);

        Assert.assertEquals(hdfsInput.getFeaturePath(), expectedFeaturePath, "wrong featurePath");
        Assert.assertEquals(hdfsInput.getName(), expectedLogicalName, "wrong logical name");
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
    public void testFeatureCodecCache() {
        Assert.assertEquals(getVariantFeatureInputWithCachedCodec().getFeatureCodecClass(), VCFCodec.class);
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testFeatureCodecCacheSerialization() throws IOException, ClassNotFoundException {
        FeatureInput<VariantContext> featureInput = getVariantFeatureInputWithCachedCodec();

        // serialize
        byte[] serializedFeatureInput;
        try (final ByteArrayOutputStream bos = new ByteArrayOutputStream();
             final ObjectOutputStream oos = new ObjectOutputStream(bos)) {
            oos.writeObject(featureInput);
            serializedFeatureInput = bos.toByteArray();
        }
        Assert.assertNotNull(serializedFeatureInput);

        // deserialize
        FeatureInput<VariantContext> roundTrippedFeatureInput;
        try (final ObjectInputStream ois = new ObjectInputStream(new ByteArrayInputStream(serializedFeatureInput))) {
            roundTrippedFeatureInput = (FeatureInput<VariantContext>) ois.readObject();
        }
        Assert.assertNotNull(roundTrippedFeatureInput);

        // we expect to lose the cached feature codec class on serialization, but retain the feature path
        Assert.assertNull(roundTrippedFeatureInput.getFeatureCodecClass());
        Assert.assertEquals(featureInput.getFeaturePath(), roundTrippedFeatureInput.getFeaturePath());
    }

    @SuppressWarnings("unchecked")
    private FeatureInput<VariantContext> getVariantFeatureInputWithCachedCodec() {
        final File inputVCFFile = new File(FEATURE_INPUT_TEST_DIRECTORY, "minimal_vcf4_file.vcf");
        final FeatureInput<VariantContext> featureInput = new FeatureInput<>(inputVCFFile.getAbsolutePath());
        Assert.assertNull(featureInput.getFeatureCodecClass());

        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(new File(featureInput.getFeaturePath()));
        featureInput.setFeatureCodecClass((Class<FeatureCodec<VariantContext, ?>>)codec.getClass());

        return featureInput;
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

    @DataProvider(name = "HasUserSuppliedNameData")
    public Object[][] hasUserSuppliedNameData() {
        return new Object[][] {
                {"hdfs://localhost/user/my.vcf", false},
                {"myname:hdfs://localhost/user/my.vcf", true},
                {"myname,key1=value1:hdfs://localhost/user/my.vcf", true},
                {"myname//:hdfs://localhost/user/my.vcf", true},
                {"myname//:/user/my.vcf", true},
                {"/user/my.vcf", false},
        };
    }

    @Test(dataProvider = "HasUserSuppliedNameData")
    public void testHasUserSuppliedName(final String inputString, final boolean isUserSupplied) {
        final FeatureInput<VariantContext> input = new FeatureInput<>(inputString);
        Assert.assertEquals(input.hasUserSuppliedName(), isUserSupplied);
    }

}
