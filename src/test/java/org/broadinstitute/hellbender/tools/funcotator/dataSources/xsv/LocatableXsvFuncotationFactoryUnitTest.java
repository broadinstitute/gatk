package org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv;

import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.FeatureInputTestTools;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationBuilder;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvLocatableTableCodec;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

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
                                                          final List<Funcotation> expected,
                                                          final boolean doGzip) {
        return new Object[] {
            createVariantContext(contig, start, end, refAlleleString, altAlleleString, refFilePath),
            new ReferenceContext( referenceDataSourceMap.get(refFilePath), new SimpleInterval(contig, start, end)),
            funcotationFieldNames,
            featureList,
            gencodeFuncotationList,
            expected,
            doGzip
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
        Object[][] arr1 = provideDataForTestCreateFuncotations(false);
        Object[][] arr2 = provideDataForTestCreateFuncotations(true);
        return Utils.concat(arr1, arr2, Object[][]::new);
    }

    private Object[][] provideDataForTestCreateFuncotations(boolean doGzip) {

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
                        Collections.singletonList(TableFuncotation.create(reportableFieldNames, emptyFieldList, defaultAltAllele, defaultDataSourceName, null)),
                        doGzip
                ),
                // Trivial case where null Features are in the list:
                helpProvideForTestCreateFuncotations(
                        "chr3", 178866314, 178866314,
                        "C", defaultAltAllele.getBaseString(), FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref(),
                        reportableFieldNames,
                        Arrays.asList(null, null, null), Collections.emptyList(),
                        Collections.singletonList(TableFuncotation.create(reportableFieldNames, emptyFieldList, defaultAltAllele, defaultDataSourceName, null)),
                        doGzip
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
                        Collections.singletonList(TableFuncotation.create(xsvTableFeature1, defaultAltAllele, defaultDataSourceName, null)),
                        doGzip
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
                        Collections.singletonList(TableFuncotation.create(xsvTableFeature1, defaultAltAllele, defaultDataSourceName, null)),
                        doGzip
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
                        Collections.singletonList(TableFuncotation.create(xsvTableFeature1, defaultAltAllele, defaultDataSourceName, null)),
                        doGzip
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
                        Collections.singletonList(TableFuncotation.create(xsvTableFeature1, defaultAltAllele, defaultDataSourceName, null)),
                        doGzip
                ),
        };
    }

    @DataProvider
    private Object[][] provideForTestSetSupportedFuncotationFields() {
        return new Object[][] {
                // One Valid XSV (csv) Locatable Data File:
                {
                        IOUtils.getPath(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE1_DATA_PATH),
                        IOUtils.getPath(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE1_CONFIG_PATH),
                        new LinkedHashSet<>(Arrays.asList("XSV_LOCATABLE_TEST_NAME_Villain", "XSV_LOCATABLE_TEST_NAME_test_val", "XSV_LOCATABLE_TEST_NAME_Bond"))
                },
                // One Valid XSV (csv) Locatable Data File:
                {
                        IOUtils.getPath(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE2_DATA_PATH),
                        IOUtils.getPath(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE2_CONFIG_PATH),
                        new LinkedHashSet<>(Arrays.asList("SECOND_XSV_NAME_Car_Maker", "SECOND_XSV_NAME_Tire_Maker", "SECOND_XSV_NAME_Parent_Company"))
                },
                // One Valid XSV (tsv) Locatable Data File:
                {
                        IOUtils.getPath(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE3_DATA_PATH),
                        IOUtils.getPath(FuncotatorTestConstants.XSV_LOCATABLE_TEST_FILE3_CONFIG_PATH),
                        new LinkedHashSet<>(Arrays.asList("XSV_LOCATABLE_TEST_NAME_Villain", "XSV_LOCATABLE_TEST_NAME_test_val", "XSV_LOCATABLE_TEST_NAME_Bond"))
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test
    public void testRequiresFeatures() {
        Assert.assertTrue(new LocatableXsvFuncotationFactory(LocatableXsvFuncotationFactory.DEFAULT_NAME, DataSourceFuncotationFactory.DEFAULT_VERSION_STRING, new LinkedHashMap<>(), null).requiresFeatures());
    }

    @Test(dataProvider = "provideForTestGetName")
    public void testGetName(final String name, final String expected) {
        final LocatableXsvFuncotationFactory locatableXsvFuncotationFactory;
        if ( name == null ) {
            locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory(LocatableXsvFuncotationFactory.DEFAULT_NAME, DataSourceFuncotationFactory.DEFAULT_VERSION_STRING, new LinkedHashMap<>(), null);
        }
        else {
            locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory(name, DataSourceFuncotationFactory.DEFAULT_VERSION_STRING, new LinkedHashMap<>(), null);
        }

        Assert.assertEquals( locatableXsvFuncotationFactory.getName(), expected );
    }

    @Test(dataProvider = "provideForTestCreateFuncotations")
    public void testCreateFuncotations(final VariantContext variant,
                                       final ReferenceContext referenceContext,
                                       final List<String> reportableFuncotationFieldNames,
                                       final List<Feature> featureList,
                                       final List<GencodeFuncotation> gencodeFuncotations,
                                       final List<Funcotation> expected,
                                       final boolean doGzip) {

        // Create a temporary file for the "backing data" which will only contain the header:
        final Path headerBackingDataFilePath = createTempPath("headerBackingDataFile", "csv" + (doGzip ? ".gz" : ""));
        final Path configFilePath;
        try (BufferedWriter writer = IOUtil.openFileForBufferedUtf8Writing(headerBackingDataFilePath)){
            writer.write("CONTIG,START,END," + reportableFuncotationFieldNames.stream().collect(Collectors.joining(",")));

            // Create a temporary file for the config file that points to the temporary file for the backing data:
            configFilePath = createTemporaryConfigFile(headerBackingDataFilePath);
        }
        catch (final IOException ex) {
            throw new GATKException("Could not write to temp file for testing: " + headerBackingDataFilePath.toUri(), ex);
        }

        final FeatureInput<? extends Feature> featureInput                   = FeatureInputTestTools.createFeatureInput( configFilePath.toUri().toString(), defaultDataSourceName );
        final LocatableXsvFuncotationFactory  locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory(defaultDataSourceName, DataSourceFuncotationFactory.DEFAULT_VERSION_STRING, new LinkedHashMap<>(), featureInput);
        locatableXsvFuncotationFactory.setSupportedFuncotationFields(headerBackingDataFilePath);

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
    public void testSetSupportedFuncotationFields(final Path dataFilePath,
                                                  final Path configFilePath,
                                                  final LinkedHashSet<String> expected) {

        final FeatureInput<? extends Feature> featureInput                   = FeatureInputTestTools.createFeatureInput(configFilePath.toUri().toString(), defaultDataSourceName);
        final LocatableXsvFuncotationFactory  locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory(LocatableXsvFuncotationFactory.DEFAULT_NAME, DataSourceFuncotationFactory.DEFAULT_VERSION_STRING, new LinkedHashMap<>(), featureInput);

        locatableXsvFuncotationFactory.setSupportedFuncotationFields(dataFilePath);

        Assert.assertEquals(
                locatableXsvFuncotationFactory.getSupportedFuncotationFields(),
                expected
        );
    }

    @Test(expectedExceptions = GATKException.class)
    public void testGetSupportedFuncotationFields() {
        final LocatableXsvFuncotationFactory locatableXsvFuncotationFactory = new LocatableXsvFuncotationFactory(LocatableXsvFuncotationFactory.DEFAULT_NAME, DataSourceFuncotationFactory.DEFAULT_VERSION_STRING, new LinkedHashMap<>(), null);
        locatableXsvFuncotationFactory.getSupportedFuncotationFields();
    }

    private Path createTemporaryConfigFile(final Path backingDataSourcePath) throws IOException {
        return createTemporaryConfigFile(backingDataSourcePath, ",");
    }

    private Path createTemporaryConfigFile(final Path backingDataSourcePath, final String delimiter) throws IOException {
        // Config file must be next to backingDataSourcePath, and have the same base name, with the .config extension:
        final String backingDataSourceFileName = backingDataSourcePath.toFile().getName();
        final String configFileBaseName = FilenameUtils.removeExtension(backingDataSourceFileName);
        final Path configPath = backingDataSourcePath.resolveSibling(configFileBaseName + XsvLocatableTableCodec.CONFIG_FILE_EXTENSION);

        final File configFile = configPath.toAbsolutePath().toFile();
        configFile.createNewFile();

        try(final PrintWriter writer = new PrintWriter(configPath.toAbsolutePath().toFile())) {
            writer.println("name = ");
            writer.println("version = TEST");
            writer.println("src_file = " + backingDataSourceFileName);
            writer.println("origin_location = LocatableXsvFuncotationFactoryUnitTest.java");
            writer.println("preprocessing_script = ");
            writer.println("");
            writer.println("# Supported types:");
            writer.println("# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID");
            writer.println("# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location");
            writer.println("# gencode      -- Custom datasource class for GENCODE");
            writer.println("# cosmic       -- Custom datasource class for COSMIC");
            writer.println("# vcf          -- Custom datasource class for Variant Call Format (VCF) files");
            writer.println("                    type = locatableXSV");
            writer.println("");
            writer.println("# Required field for GENCODE files.");
            writer.println("# Path to the FASTA file from which to load the sequences for GENCODE transcripts:");
            writer.println("gencode_fasta_path =");
            writer.println("");
            writer.println("# Required field for simpleXSV files.");
            writer.println("# Valid values:");
            writer.println("#     GENE_NAME");
            writer.println("#     TRANSCRIPT_ID");
            writer.println("xsv_key = ");
            writer.println("");
            writer.println("# Required field for simpleXSV files.");
            writer.println("# The 0-based index of the column containing the key on which to match");
            writer.println("                    xsv_key_column =");
            writer.println("");
            writer.println("# Required field for simpleXSV AND locatableXSV files.");
            writer.println("# The delimiter by which to split the XSV file into columns.");
            writer.println("xsv_delimiter = " + delimiter);
            writer.println("");
            writer.println("# Required field for simpleXSV files.");
            writer.println("# Whether to permissively match the number of columns in the header and data rows");
            writer.println("# Valid values:");
            writer.println("#     true");
            writer.println("#     false");
            writer.println("xsv_permissive_cols = ");
            writer.println("");
            writer.println("# Required field for locatableXSV files.");
            writer.println("# The 0-based index of the column containing the contig for each row");
            writer.println("contig_column = 0 ");
            writer.println("");
            writer.println("# Required field for locatableXSV files.");
            writer.println("# The 0-based index of the column containing the start position for each row");
            writer.println("start_column = 1 ");
            writer.println("");
            writer.println("# Required field for locatableXSV files.");
            writer.println("# The 0-based index of the column containing the end position for each row");
            writer.println("end_column = 2");
        }

        return configPath;
    }

    @Test
    public void testNoSupportOfSegments() {
        final LocatableXsvFuncotationFactory factory = new LocatableXsvFuncotationFactory(LocatableXsvFuncotationFactory.DEFAULT_NAME, DataSourceFuncotationFactory.DEFAULT_VERSION_STRING, new LinkedHashMap<>(), null);

        Assert.assertFalse(factory.isSupportingSegmentFuncotation());
        Assert.assertEquals(factory.getSupportedFuncotationFieldsForSegments(), Collections.emptyList());
    }
}
