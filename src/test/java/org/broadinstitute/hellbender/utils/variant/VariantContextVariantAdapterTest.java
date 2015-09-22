package org.broadinstitute.hellbender.utils.variant;

import com.google.cloud.dataflow.sdk.Pipeline;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.coders.VariantCoder;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class VariantContextVariantAdapterTest extends BaseTest {
    // Clearly, only for testing.
    public static VariantContextVariantAdapter createVariantContextVariantAdapterForTesting(VariantContextVariantAdapter vc, UUID uuid) {
        return new VariantContextVariantAdapter(vc, uuid);
    }

    private static final String FEATURE_DATA_SOURCE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final File QUERY_TEST_VCF = new File(FEATURE_DATA_SOURCE_TEST_DIRECTORY + "feature_data_source_test.vcf");

    UUID defaultUUID() {
        return new UUID(0L, 0L);
    }

    @Test(dataProvider = "VariantDataProvider")
    public void testVariantAdapter(final List<Variant> expectedVariantList) {
        // The test suite for reading in VCF files is FeatureDataSourceUnitTest.
        try (FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, new VCFCodec())) {
            List<Variant> variantList = new ArrayList<>();
            for (VariantContext feature : featureSource) {
                VariantContextVariantAdapter va = new VariantContextVariantAdapter(feature);
                variantList.add(va);
            }
            // Now, test to see that every variant is in in the expected set.
            Assert.assertEquals(variantList.size(), expectedVariantList.size());
            for (Variant v : variantList) {
                boolean matchFound = false;
                for (Variant vv : expectedVariantList) {
                    if (VariantUtils.variantsAreEqual(v, vv)) {
                        matchFound = true;
                    }
                }
                Assert.assertTrue(matchFound, v.toString());
            }
        }
    }

    @Test(dataProvider = "FileVariantDataProvider")
    public void coderTest(final List<Variant> variantList) {
        // The simplest way to figure out if a class is coded correctly is to create a PCollection
        // of that type and see if matches the List version.
        Pipeline p = GATKTestPipeline.create();
        DataflowTestUtils.pCollectionCreateAndVerify(p, variantList, new VariantCoder());
        p.run();
    }


    @DataProvider(name = "VariantDataProvider")
    public Object[][] getVariantData() {
        // We use the defaultUUID because we want to have be able to use equals (without clearing, we'd see items
        // aren't the same because of UUIDs).
        List<Variant> variantSet = new ArrayList<>();
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 100, 100), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 199, 200), false, true, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 200, 200), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 203, 206), false, true, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 280, 280), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 284, 286), false, true, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 285, 285), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 286, 286), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 999, 999), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 1000, 1000), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 1000, 1003), false, true, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 1076, 1076), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 1150, 1150), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 1176, 1176), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("2", 200, 200), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("2", 525, 525), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("2", 548, 550), false, true, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("2", 640, 640), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("2", 700, 700), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("3", 1, 1), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("3", 300, 300), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("3", 300, 303), false, true, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("3", 400, 400), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("4", 600, 600), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("4", 775, 775), true, false, defaultUUID()));
        variantSet.add(new MinimalVariant(new SimpleInterval("4", 776, 779), false, true, defaultUUID()));
        return new Object[][]{
                { variantSet }
        };
    }

    @DataProvider(name = "FileVariantDataProvider")
    public Object[][] getFileVariantData() {
        List<Variant> variantList = new ArrayList<>();
        try (FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, new VCFCodec())) {
            for (VariantContext feature : featureSource) {
                VariantContextVariantAdapter va = new VariantContextVariantAdapter(feature);
                variantList.add(va);
            }
        }
        return new Object[][]{
                {variantList},
        };
    }
}
