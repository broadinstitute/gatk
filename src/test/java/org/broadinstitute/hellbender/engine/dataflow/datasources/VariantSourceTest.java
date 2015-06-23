package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.SkeletonVariant;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.broadinstitute.hellbender.utils.variant.VariantUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.UUID;

public class VariantSourceTest extends BaseTest {
    private static final String FEATURE_DATA_SOURCE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final File QUERY_TEST_VCF = new File(FEATURE_DATA_SOURCE_TEST_DIRECTORY + "feature_data_source_test.vcf");

    UUID defaultUUID() {
        return new UUID(0L, 0L);
    }

    static class clearUUIDDoFn extends DoFn<Variant, Variant> {
        private static final long serialVersionUID = 1L;
        @Override
        public void processElement(ProcessContext c) throws Exception {
            ((VariantContextVariantAdapter) c.element()).clearUUID();
            c.output(c.element());
        }
    }

    @Test(dataProvider = "RealVariantData")
    public void testVariantPCollection(final List<Variant> variants) {
        List<Variant> variants2 = Lists.newArrayList(variants.iterator());

        // Now make a PCollection, to verify that variants can be coded.
        Pipeline p = TestPipeline.create();

        PCollection<Variant> pVariants = p.apply(Create.of(variants));
        DataflowAssert.that(pVariants).containsInAnyOrder(variants2);

        // Now, we can test that we can get a PCollection from a file
        List<String> files = Arrays.asList(FEATURE_DATA_SOURCE_TEST_DIRECTORY + "feature_data_source_test.vcf");
        VariantsDataflowSource variantsDataflowSource = new VariantsDataflowSource(files, p);
        PCollection<Variant> allVariants = variantsDataflowSource.getAllVariants();


        // We have to clear the UUIDs to make the comparison.
        allVariants = allVariants.apply(ParDo.of(new clearUUIDDoFn()));
        for (Variant v : variants2) {
            ((VariantContextVariantAdapter) v).clearUUID();
        }
        DataflowAssert.that(allVariants).containsInAnyOrder(variants2);

        p.run();
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
        Pipeline p = TestPipeline.create();
        // Turns out we don't need to register a coder for VariantContextVariantAdapter.

        PCollection<Variant> pShards = p.apply(Create.of(variantList));

        List<Variant> sameVariants = Lists.newArrayList();
        Assert.assertTrue(sameVariants.addAll(variantList));
        DataflowAssert.that(pShards).containsInAnyOrder(sameVariants);
        p.run();
    }


    @DataProvider(name = "VariantDataProvider")
    public Object[][] getVariantData() {
        List<Variant> variantSet = new ArrayList<>();
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 100, 100), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 199, 200), false, true, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 200, 200), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 203, 206), false, true, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 280, 280), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 284, 286), false, true, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 285, 285), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 286, 286), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 999, 999), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 1000, 1000), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 1000, 1003), false, true, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 1076, 1076), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 1150, 1150), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("1", 1176, 1176), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("2", 200, 200), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("2", 525, 525), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("2", 548, 550), false, true, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("2", 640, 640), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("2", 700, 700), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("3", 1, 1), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("3", 300, 300), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("3", 300, 303), false, true, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("3", 400, 400), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("4", 600, 600), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("4", 775, 775), true, false, defaultUUID()));
        variantSet.add(new SkeletonVariant(new SimpleInterval("4", 776, 779), false, true, defaultUUID()));
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
    @DataProvider(name = "RealVariantData")
    public Object[][] getRealVariantData() {
        List<String> files = Lists.newArrayList(FEATURE_DATA_SOURCE_TEST_DIRECTORY + "feature_data_source_test.vcf");
        List<Variant> variants = VariantsDataflowSource.getVariantsList(files);
        return new Object[][]{
            {variants},
        };

    }
}
