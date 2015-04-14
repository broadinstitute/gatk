package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCF3Codec;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.table.TableCodec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class FeatureManagerUnitTest extends BaseTest {
    private static final String FEATURE_MANAGER_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

    @DataProvider(name = "DetectCorrectFileFormatTestData")
    public Object[][] getDetectCorrectFileFormatTestData() {
        return new Object[][] {
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_vcf4_file.vcf"), VCFCodec.class },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_vcf3_file.vcf"), VCF3Codec.class },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_bcf_file.bcf"), BCF2Codec.class },
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_bed_file.bed"), BEDCodec.class},
                { new File(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_table_file.table"), TableCodec.class}
        };
    }

    @Test(dataProvider = "DetectCorrectFileFormatTestData")
    public void testDetectCorrectFileFormat( final File file, final Class<?> expectedCodecClass ) {
        Assert.assertEquals(FeatureManager.getCodecForFile(file).getClass(), expectedCodecClass,
                            "Wrong codec selected for file " + file.getAbsolutePath());
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testDetectUnsupportedFileFormat() {
        final File unsupportedFile = new File(FEATURE_MANAGER_TEST_DIRECTORY + "unsupported_format_file");

        // Eliminate file existence/readability as a factor in the test
        Assert.assertTrue(unsupportedFile.canRead(), "Cannot test detection of unsupported file formats on an unreadable file");

        // Should throw, since the file exists and is readable, but is in an unsupported format
        FeatureManager.getCodecForFile(unsupportedFile);
    }

    @CommandLineProgramProperties(usage = "", usageShort = "")
    private static class ValidFeatureArgumentSource extends CommandLineProgram {
        // We should be able to detect the type parameter of a FeatureInput regardless of whether or
        // not it's wrapped within a Collection

        @Argument(fullName = "variantContextFeatureInput")
        public FeatureInput<VariantContext> variantContextFeatureInput;

        @Argument(fullName = "variantContextListFeatureInput")
        public List<FeatureInput<VariantContext>> variantContextListFeatureInput;

        @Argument(fullName = "bedFeatureInput")
        public FeatureInput<BEDFeature> bedFeatureInput;

        @Argument(fullName = "bedListFeatureInput")
        public List<FeatureInput<BEDFeature>> bedListFeatureInput;

        public ValidFeatureArgumentSource() {
            // The real argument parsing system ensures that Collection argument fields are non-null
            variantContextListFeatureInput = new ArrayList<>();
            bedListFeatureInput = new ArrayList<>();
        }

        @Override
        protected Object doWork() {
            return null;
        }
    }

    @CommandLineProgramProperties(usage = "", usageShort = "")
    @SuppressWarnings("rawtypes")
    private static class InvalidFeatureArgumentSource extends CommandLineProgram {
        // FeatureInputs without type parameters (ie., raw types) should be detected as errors

        @Argument(fullName = "parameterlessFeatureInput")
        public FeatureInput parameterlessFeatureInput;

        @Argument(fullName = "parameterlessListFeatureInput")
        public List<FeatureInput> parameterlessListFeatureInput;

        @Override
        protected Object doWork() {
            return null;
        }
    }

    @DataProvider(name = "DetectFeatureTypeFromFeatureInputFieldTestData")
    public Object[][] getDetectFeatureTypeFromFeatureInputFieldTestData() {
        return new Object[][] {
                { "variantContextFeatureInput", VariantContext.class },
                { "variantContextListFeatureInput", VariantContext.class },
                { "bedFeatureInput", BEDFeature.class },
                { "bedListFeatureInput", BEDFeature.class }
        };
    }

    @Test(dataProvider = "DetectFeatureTypeFromFeatureInputFieldTestData")
    public void testDetectFeatureTypeFromFeatureInputField( final String fieldName, final Class<?> expectedFeatureType ) throws NoSuchFieldException {
        Assert.assertEquals(FeatureManager.getFeatureTypeForFeatureInputField(ValidFeatureArgumentSource.class.getDeclaredField(fieldName)),
                            expectedFeatureType,
                            "Wrong Feature type detected for field " + fieldName);
    }

    @DataProvider(name = "DetectParameterlessFeatureInputsTestData")
    public Object[][] getDetectParameterlessFeatureInputsTestData() {
        return new Object[][] {
                { "parameterlessFeatureInput" },
                { "parameterlessListFeatureInput" }
        };
    }

    @Test(dataProvider = "DetectParameterlessFeatureInputsTestData", expectedExceptions = GATKException.class)
    public void testDetectParameterlessFeatureInputs( final String fieldName ) throws NoSuchFieldException {
        // FeatureInput fields that lack a type parameter should cause a GATKException to be thrown
        FeatureManager.getFeatureTypeForFeatureInputField(InvalidFeatureArgumentSource.class.getDeclaredField(fieldName));
    }

    @Test
    public void testHandleRequestForValidFeatureInputs() {
        ValidFeatureArgumentSource toolInstance = new ValidFeatureArgumentSource();

        // Initialize two of the FeatureInput fields as they would be initialized by the argument-parsing
        // system to simulate a run of the tool with two FeatureInputs.
        toolInstance.variantContextFeatureInput = new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "feature_data_source_test.vcf");
        toolInstance.bedListFeatureInput.add(new FeatureInput<>(FEATURE_MANAGER_TEST_DIRECTORY + "minimal_bed_file.bed"));

        FeatureManager manager = new FeatureManager(toolInstance);
        List<VariantContext> vcFeatures = manager.getFeatures(toolInstance.variantContextFeatureInput, new SimpleInterval("1", 1, 2000));
        List<BEDFeature> bedFeatures = manager.getFeatures(toolInstance.bedListFeatureInput.get(0), new SimpleInterval("1", 1, 1));

        Assert.assertEquals(vcFeatures.size(), 14, "Wrong number of Features returned from VariantContext test Feature file");
        Assert.assertEquals(bedFeatures.size(), 1, "Wrong number of Features returned from BED test Feature file");
    }

    @Test(expectedExceptions = GATKException.class)
    public void testHandleRequestForNonExistentFeatureInput() {
        FeatureManager manager = new FeatureManager(new ValidFeatureArgumentSource());

        // Requests for FeatureInputs not declared in the tool's class hierarchy (or associated ArgumentCollections)
        // should throw an exception
        manager.getFeatures(new FeatureInput<Feature>("featureInputNotDeclaredInTool"), new SimpleInterval("1", 1, 1));
    }
}
