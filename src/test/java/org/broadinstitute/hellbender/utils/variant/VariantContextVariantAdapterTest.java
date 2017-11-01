package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class VariantContextVariantAdapterTest extends GATKBaseTest {

    private static final String FEATURE_DATA_SOURCE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final File QUERY_TEST_VCF = new File(FEATURE_DATA_SOURCE_TEST_DIRECTORY + "feature_data_source_test.vcf");

    @Test(dataProvider = "VariantDataProvider")
    public void testVariantAdapter(final List<GATKVariant> expectedVariantList) {
        // The test suite for reading in VCF files is FeatureDataSourceUnitTest.
        try (FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF)) {
            List<GATKVariant> variantList = new ArrayList<>();
            for (VariantContext feature : featureSource) {
                VariantContextVariantAdapter va = new VariantContextVariantAdapter(feature);
                variantList.add(va);
            }
            // Now, test to see that every variant is in in the expected set.
            Assert.assertEquals(variantList.size(), expectedVariantList.size());
            for (GATKVariant v : variantList) {
                boolean matchFound = false;
                for (GATKVariant vv : expectedVariantList) {
                    if (VariantUtils.variantsAreEqual(v, vv)) {
                        matchFound = true;
                    }
                }
                Assert.assertTrue(matchFound, v.toString());
            }
        }
    }

    @DataProvider(name = "VariantDataProvider")
    public Object[][] getVariantData() {
        List<GATKVariant> variantSet = new ArrayList<>();
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 100, 100), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 199, 200), false, true));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 200, 200), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 203, 206), false, true));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 280, 280), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 284, 286), false, true));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 285, 285), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 286, 286), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 999, 999), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 1000, 1000), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 1000, 1003), false, true));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 1076, 1076), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 1150, 1150), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("1", 1176, 1176), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("2", 200, 200), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("2", 525, 525), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("2", 548, 550), false, true));
        variantSet.add(new MinimalVariant(new SimpleInterval("2", 640, 640), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("2", 700, 700), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("3", 1, 1), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("3", 300, 300), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("3", 300, 303), false, true));
        variantSet.add(new MinimalVariant(new SimpleInterval("3", 400, 400), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("4", 600, 600), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("4", 775, 775), true, false));
        variantSet.add(new MinimalVariant(new SimpleInterval("4", 776, 779), false, true));
        return new Object[][]{
                { variantSet }
        };
    }
}
