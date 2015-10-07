package org.broadinstitute.hellbender.engine.dataflow.coders;

import com.google.cloud.dataflow.sdk.Pipeline;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class VariantContextVariantAdapterCoderTest extends BaseTest {

    private static final String FEATURE_DATA_SOURCE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final File QUERY_TEST_VCF = new File(FEATURE_DATA_SOURCE_TEST_DIRECTORY + "feature_data_source_test.vcf");

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

    @Test(dataProvider = "FileVariantDataProvider")
    public void coderTest(final List<Variant> variantList) {
        // The simplest way to figure out if a class is coded correctly is to create a PCollection
        // of that type and see if matches the List version.
        Pipeline p = GATKTestPipeline.create();
        DataflowTestUtils.pCollectionCreateAndVerify(p, variantList, new VariantCoder());
        p.run();
    }
}
