package org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationBuilder;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

/**
 * Unit tests for {@link LocatableXsvFuncotationFactory}.
 * Created by jonn on 12/11/17.
 */
public class LocatableXsvFuncotationFactoryUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    private static final Allele defaultAltAllele = Allele.create("G", false);
    private static final String defaultDataSourceName = "LocatableXsvFuncotationFactoryUnitTest";
    private static final Map<String, ReferenceDataSource> referenceDataSourceMap;

    static {
        referenceDataSourceMap = new HashMap<>(2);

        referenceDataSourceMap.put(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(), ReferenceDataSource.of( new File(FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref()).toPath() ));
        referenceDataSourceMap.put(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(), ReferenceDataSource.of( new File(FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()).toPath() ));
    }

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Data Types:
    private static class DummyTestFeature implements Feature {

        private final String contig;
        private final int start;
        private final int stop;

        public DummyTestFeature(final String contig, final int start, final int stop) {
            this.contig = contig;
            this.start = start;
            this.stop = stop;
        }

        @Override
        public String getContig() {
            return contig;
        }

        @Override
        public int getStart() {
            return start;
        }

        @Override
        public int getEnd() {
            return stop;
        }
    }

    //==================================================================================================================
    // Helper Methods:

    private VariantContext createVariantContext(final String contig,
                                                final int start,
                                                final int end,
                                                final String refString,
                                                final String altString,
                                                final String referenceFilePath) {

        final Allele refAllele = Allele.create(refString, true);
        final Allele altAllele = Allele.create(altString);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                referenceFilePath,
                contig,
                start,
                end,
                Arrays.asList(refAllele, altAllele)
        );

        return variantContextBuilder.make();
    }

    private Object[] helpProvideForTestCreateFuncotations(final String contig,
                                                          final int start,
                                                          final int end,
                                                          final String refAlleleString,
                                                          final String altAlleleString,
                                                          final String refFilePath,
                                                          final List<String> funcotationFieldNames,
                                                          final List<Feature> featureList,
                                                          final List<GencodeFuncotation> gencodeFuncotationList,
                                                          final List<Funcotation> expected) {
        return new Object[] {
            createVariantContext(contig, start, end, refAlleleString, altAlleleString, refFilePath),
            new ReferenceContext( referenceDataSourceMap.get(refFilePath), new SimpleInterval(contig, start, end)),
            funcotationFieldNames,
            featureList,
            gencodeFuncotationList,
            expected
        };
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Object[][] provideForTestGetName() {
        return new Object[][] {
                // No name provided / default name:
                {null, "LocatableXsv"},
                // Any name string:
                {"TEST", "TEST"},
                // Any name string:
                {"4", "4"},
        };
    }

    @DataProvider
    private Object[][] provideForTestCreateFuncotations() {

        final List<String> fieldNames = Arrays.asList("0", "1", "2", "3", "4", "5", "6");

        // Fields that contain data that will be reported (fields that do not notionally contain location info):
        final List<String> reportableFieldNames = fieldNames.subList(3,7);
        final List<String> emptyFieldList = Arrays.asList("", "", "", "");

        final XsvTableFeature xsvTableFeature1 = new XsvTableFeature(
                0,1,2,
                fieldNames,
                Arrays.asList("chr3", "178866310", "178866320", "this", "is", "a", "test"),
                defaultDataSourceName
        );

        final XsvTableFeature xsvTableFeature2 = new XsvTableFeature(
                0,1,2,
                fieldNames,
                Arrays.asList("chr3", "178866300", "178866350", "this", "is", "test", "2"),
                defaultDataSourceName
        );

        final XsvTableFeature xsvTableFeature3 = new XsvTableFeature(
                0,1,2,
                fieldNames,
                Arrays.asList("chr3", "178866000", "178866550", "this", "is", "test", "3"),
                defaultDataSourceName
        );

        return new Object[][] {
                // Trivial case the list of features is empty:
                helpProvideForTestCreateFuncotations(
                        "chr3", 178866314, 178866314,
                        "C", defaultAltAllele.getBaseString(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        reportableFieldNames,
                        Collections.emptyList(), Collections.emptyList(),
                        Collections.singletonList(TableFuncotation.create(reportableFieldNames, emptyFieldList, defaultAltAllele, defaultDataSourceName, null))
                ),
                // Trivial case where null Features are in the list:
                helpProvideForTestCreateFuncotations(
                        "chr3", 178866314, 178866314,
                        "C", defaultAltAllele.getBaseString(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        reportableFieldNames,
                        Arrays.asList(null, null, null), Collections.emptyList(),
                        Collections.singletonList(TableFuncotation.create(reportableFieldNames, emptyFieldList, defaultAltAllele, defaultDataSourceName, null))
                ),
                // One XsvTableFeature in list
                helpProvideForTestCreateFuncotations(
                        "chr3", 178866314, 178866314,
                        "C", defaultAltAllele.getBaseString(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        fieldNames,
                        Collections.singletonList(
                            xsvTableFeature1
                        ),
                        Collections.emptyList(),
                        Collections.singletonList(TableFuncotation.create(xsvTableFeature1, defaultAltAllele, defaultDataSourceName, null))
                ),
                // Two XsvTableFeatures in list
                helpProvideForTestCreateFuncotations(
                        "chr3", 178866314, 178866314,
                        "C", defaultAltAllele.getBaseString(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        fieldNames,
                        Arrays.asList(
                                xsvTableFeature1, xsvTableFeature2
                        ),
                        Collections.emptyList(),
                        // TODO: Commented out because of issue #4930.  When issue is fixed, revert to this test case! (https://github.com/broadinstitute/gatk/issues/4930)
//                        Arrays.asList(TableFuncotation.create(xsvTableFeature1, defaultAltAllele, defaultDataSourceName, null), TableFuncotation.create(xsvTableFeature2, defaultAltAllele, defaultDataSourceName, null))
                        Collections.singletonList(TableFuncotation.create(xsvTableFeature1, defaultAltAllele, defaultDataSourceName, null))
                ),
                // Many XsvTableFeatures in list
                helpProvideForTestCreateFuncotations(
                        "chr3", 178866314, 178866314,
                        "C", defaultAltAllele.getBaseString(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        fieldNames,
                        Arrays.asList(
                                xsvTableFeature1, xsvTableFeature2, xsvTableFeature3
                        ),
                        Collections.emptyList(),
                        // TODO: Commented out because of issue #4930.  When issue is fixed, revert to this test case! (https://github.com/broadinstitute/gatk/issues/4930)
//                        Arrays.asList(TableFuncotation.create(xsvTableFeature1, defaultAltAllele, defaultDataSourceName, null), TableFuncotation.create(xsvTableFeature2, defaultAltAllele, defaultDataSourceName, null), TableFuncotation.create(xsvTableFeature3, defaultAltAllele, defaultDataSourceName, null))
                        Collections.singletonList(TableFuncotation.create(xsvTableFeature1, defaultAltAllele, defaultDataSourceName, null))
                ),
                // Many XsvTableFeatures in list and non-empty GencodeFuncotations
                helpProvideForTestCreateFuncotations(
                        "chr3", 178866314, 178866314,
                        "C", defaultAltAllele.getBaseString(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        fieldNames,
                        Arrays.asList(
                                xsvTableFeature1, xsvTableFeature2, xsvTableFeature3
                        ),
                        Collections.singletonList(
                                new GencodeFuncotationBuilder().setChromosome("chr3").setStart(178866314).setEnd(178866314).build()
                        ),
                        // TODO: Commented out because of issue #4930.  When issue is fixed, revert to this test case! (https://github.com/broadinstitute/gatk/issues/4930)
//                        Arrays.asList(TableFuncotation.create(xsvTableFeature1, defaultAltAllele, defaultDataSourceName, null), TableFuncotation.create(xsvTableFeature2, defaultAltAllele, defaultDataSourceName, null), TableFuncotation.create(xsvTableFeature3, defaultAltAllele, defaultDataSourceName, null))
                        Collections.singletonList(TableFuncotation.create(xsvTableFeature1, defaultAltAllele, defaultDataSourceName, null))
                ),
        };
    }

    @DataProvider
    private Object[][] provideForTestSetSupportedFuncotationFields() {
        return new Object[][] {
                // Empty list of data files:
                {Collections.emptyList(), new LinkedHashSet<>()},
                // One Valid XSV (csv) Locatable Data File:
                {
                        Collections.singletonList(Paths.get(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE1_PATH)),
                        new LinkedHashSet<>(Arrays.asList("XSV_LOCATABLE_TEST_NAME_Villain", "XSV_LOCATABLE_TEST_NAME_test_val", "XSV_LOCATABLE_TEST_NAME_Bond"))
                },
                // One Valid XSV (tsv) Locatable Data File:
                {
                        Collections.singletonList(Paths.get(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE3_PATH)),
                        new LinkedHashSet<>(Arrays.asList("XSV_LOCATABLE_TEST_NAME_Villain", "XSV_LOCATABLE_TEST_NAME_test_val", "XSV_LOCATABLE_TEST_NAME_Bond"))
                },
                // Two Valid XSV Locatable Data Files:
                {
                    Arrays.asList(Paths.get(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE1_PATH), Paths.get(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE2_PATH)),
                    new LinkedHashSet<>(Arrays.asList("XSV_LOCATABLE_TEST_NAME_Villain", "XSV_LOCATABLE_TEST_NAME_test_val", "XSV_LOCATABLE_TEST_NAME_Bond",
                            "SECOND_XSV_NAME_Car_Maker", "SECOND_XSV_NAME_Tire_Maker", "SECOND_XSV_NAME_Parent_Company"))
                },
                // Three Valid XSV Locatable Data Files:
                {
                        Arrays.asList(Paths.get(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE1_PATH), Paths.get(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE2_PATH), Paths.get(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE3_PATH)),
                        new LinkedHashSet<>(Arrays.asList("XSV_LOCATABLE_TEST_NAME_Villain", "XSV_LOCATABLE_TEST_NAME_test_val", "XSV_LOCATABLE_TEST_NAME_Bond",
                                "SECOND_XSV_NAME_Car_Maker", "SECOND_XSV_NAME_Tire_Maker", "SECOND_XSV_NAME_Parent_Company"))
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestGetName")
    public void testGetName(final String name, final String expected) {
        final LocatableXsvFuncotationFactory locatableXsvFuncotationFactory;
        if ( name == null ) {
            locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory(LocatableXsvFuncotationFactory.DEFAULT_NAME, DataSourceFuncotationFactory.DEFAULT_VERSION_STRING);
        }
        else {
            locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory(name, DataSourceFuncotationFactory.DEFAULT_VERSION_STRING);
        }

        Assert.assertEquals( locatableXsvFuncotationFactory.getName(), expected );
    }

    @Test(dataProvider = "provideForTestCreateFuncotations")
    public void testCreateFuncotations(final VariantContext variant,
                                       final ReferenceContext referenceContext,
                                       final List<String> reportableFuncotationFieldNames,
                                       final List<Feature> featureList,
                                       final List<GencodeFuncotation> gencodeFuncotations,
                                       final List<Funcotation> expected) {

        final LocatableXsvFuncotationFactory locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory(defaultDataSourceName, DataSourceFuncotationFactory.DEFAULT_VERSION_STRING, reportableFuncotationFieldNames);


        Assert.assertEquals(
                locatableXsvFuncotationFactory.createFuncotationsOnVariant( variant, referenceContext, featureList ),
                expected
        );

        Assert.assertEquals(
                locatableXsvFuncotationFactory.createFuncotationsOnVariant( variant, referenceContext, featureList, gencodeFuncotations ),
                expected
        );
    }

    @Test(dataProvider = "provideForTestSetSupportedFuncotationFields")
    public void testSetSupportedFuncotationFields(final List<Path> dataFilePaths,
                                                  final LinkedHashSet<String> expected) {
        final LocatableXsvFuncotationFactory locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory();

        locatableXsvFuncotationFactory.setSupportedFuncotationFields(dataFilePaths);

        Assert.assertEquals(
                locatableXsvFuncotationFactory.getSupportedFuncotationFields(),
                expected
        );
    }

    @Test(expectedExceptions = GATKException.class)
    public void testGetSupportedFuncotationFields() {
        final LocatableXsvFuncotationFactory locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory();
        locatableXsvFuncotationFactory.getSupportedFuncotationFields();
    }
}
