package org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationBuilder;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.testutils.FuncotatorReferenceTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Path;
import java.util.*;

/**
 * Unit test class for {@link SimpleKeyXsvFuncotationFactory}.
 * Created by jonn on 11/28/17.
 */
public class SimpleKeyXsvFuncotationFactoryUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    private static final int squareSize = 20;
    private static final String defaultName = "XSVCSV";
    private static final List<List<String>> headerRowTable;
    private static final List<List<String>> dataTable;

    // Default locus is in PIK3CA:
    private static final String defaultContig = "chr3";
    private static final int defaultStart = 178921337;
    private static final int defaultEnd = 178921337;
    private static final Allele defaultRefAllele = Allele.create("A", true);
    private static final Allele defaultAltAllele = Allele.create("T");

    private static final VariantContext defaultVariantContext;
    private static final ReferenceContext defaultReferenceContext;

    static {
        // Initialize Static members:
        headerRowTable = new ArrayList<>(squareSize);
        dataTable = new ArrayList<>(squareSize);
        for ( int i = 0; i < squareSize; ++i ) {
            final List<String> headerRow = new ArrayList<>(squareSize);
            final List<String> dataRow = new ArrayList<>(squareSize);
            for ( int j = 0; j < squareSize; ++j ) {
                headerRow.add( defaultName + "_R" + (i+1) + "C"+ (j+1) );
                dataRow.add( "R" + (i+1) + "C"+ (j+1) );
            }
            headerRowTable.add(headerRow);
            dataTable.add(dataRow);
        }

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(
                FuncotatorReferenceTestUtils.retrieveHg19Chr19Ref(),
                defaultContig,
                defaultStart,
                defaultEnd,
                Arrays.asList(defaultRefAllele, defaultAltAllele)
        );
        defaultVariantContext = variantContextBuilder.make();

        defaultReferenceContext = new ReferenceContext(
                ReferenceDataSource.of( new File (FuncotatorReferenceTestUtils.retrieveHg19Chr3Ref()).toPath() ),
                new SimpleInterval(defaultContig, defaultStart, defaultEnd)
        );
    }

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    private <T> List<T> removeHelper(final List<T> list, final int index) {
        final ArrayList<T> tmpList = new ArrayList<>(list);

        tmpList.remove(index);

        return tmpList;
    }

    private void helpPopulateForGetSupportedFuncotationFields(final List<Object[]> outList,
                                                              final List<Integer> keyColumns,
                                                              final String path,
                                                              final String delim) {
        for ( int i = 0; i < squareSize; ++i ) {
            for ( final int keyCol : keyColumns ) {

                final List<String> headers = new ArrayList<>(headerRowTable.get(i));
                headers.remove(keyCol);

                outList.add(
                        new Object[]{
                                IOUtils.getPath(path), i, defaultName, delim, keyCol, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME,
                                new LinkedHashSet<>(headers)
                        }
                );
            }
        }
    }

    private void helpPopulateDataForCreateFuncotations(final List<Object[]> outList,
                                                       final String path,
                                                       final String delim) {

        // Programmatically make our entries:
        for ( int keyColumn = 0 ; keyColumn < squareSize ; ++keyColumn ) {
            for ( int startingHeaderRow = 0 ; startingHeaderRow < squareSize - 1 ; ++startingHeaderRow ) {

                // Add an entry for the Gene Name Key:
                outList.add(
                        new Object[]{
                                new SimpleKeyXsvFuncotationFactory(
                                        defaultName,
                                        IOUtils.getPath(path),
                                        "VERSION",
                                        delim,
                                        keyColumn,
                                        SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME,
                                        new LinkedHashMap<>(),
                                        startingHeaderRow
                                ),
                                Collections.singletonList(
                                        new GencodeFuncotationBuilder().setHugoSymbol(dataTable.get(startingHeaderRow+1).get(keyColumn)).build()
                                ),
                                Collections.singletonList(
                                        TableFuncotation.create(
                                                removeHelper(headerRowTable.get(startingHeaderRow), keyColumn),
                                                removeHelper(dataTable.get(startingHeaderRow+1), keyColumn),
                                                defaultAltAllele,
                                                defaultName, null
                                        )
                                )
                        }
                );

                // Add an entry for the Transcript ID Key:
                outList.add(
                        new Object[]{
                                new SimpleKeyXsvFuncotationFactory(
                                        defaultName,
                                        IOUtils.getPath(path),
                                        "VERSION",
                                        delim,
                                        keyColumn,
                                        SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID,
                                        new LinkedHashMap<>(),
                                        startingHeaderRow
                                ),
                                Collections.singletonList(
                                        new GencodeFuncotationBuilder().setAnnotationTranscript(dataTable.get(startingHeaderRow+1).get(keyColumn)).build()
                                ),
                                Collections.singletonList(
                                        TableFuncotation.create(
                                                removeHelper(headerRowTable.get(startingHeaderRow), keyColumn),
                                                removeHelper(dataTable.get(startingHeaderRow+1), keyColumn),
                                                defaultAltAllele,
                                                defaultName, null
                                        )
                                )
                        }
                );
            }
        }
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Iterator<Object[]> provideForTestGetSupportedFuncotationFields() {

        final List<Object[]> outList = new ArrayList<>();

        final List<Integer> keyColumnChoices = Arrays.asList(1,5,17);

        helpPopulateForGetSupportedFuncotationFields(outList, keyColumnChoices, FuncotatorTestConstants.XSV_CSV_FILE_PATH, ",");
        helpPopulateForGetSupportedFuncotationFields(outList, keyColumnChoices, FuncotatorTestConstants.XSV_TSV_FILE_PATH, "\t");
        helpPopulateForGetSupportedFuncotationFields(outList, keyColumnChoices, FuncotatorTestConstants.XSV_DEADBEEFSV_FILE_PATH, "DEADBEEF");

        return outList.iterator();
    }

    @DataProvider
    Object[][] provideForGetName() {
        return new Object[][] {
                {
                    new SimpleKeyXsvFuncotationFactory(
                            defaultName,
                            IOUtils.getPath(FuncotatorTestConstants.XSV_CSV_FILE_PATH),
                            "VERSION",
                            ",",
                            0,
                            SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME
                    ),
                    defaultName
                },
                {
                        new SimpleKeyXsvFuncotationFactory(
                                "Donatello",
                                IOUtils.getPath(FuncotatorTestConstants.XSV_CSV_FILE_PATH),
                                "VERSION",
                                ",",
                                0,
                                SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME
                        ),
                        "Donatello"
                },
                {
                        new SimpleKeyXsvFuncotationFactory(
                                "Leonardo",
                                IOUtils.getPath(FuncotatorTestConstants.XSV_CSV_FILE_PATH),
                                "VERSION",
                                ",",
                                0,
                                SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME
                        ),
                        "Leonardo"
                },
                {
                        new SimpleKeyXsvFuncotationFactory(
                                "Michelangelo",
                                IOUtils.getPath(FuncotatorTestConstants.XSV_CSV_FILE_PATH),
                                "VERSION",
                                ",",
                                0,
                                SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME
                        ),
                        "Michelangelo"
                },
                {
                new SimpleKeyXsvFuncotationFactory(
                        "Raphael",
                        IOUtils.getPath(FuncotatorTestConstants.XSV_CSV_FILE_PATH),
                        "VERSION",
                        ",",
                        0,
                        SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME
                ),
                "Raphael"
                },
        };
    }

    @DataProvider
    Iterator<Object[]> provideDataForCreateFuncotations() {

        final List<Object[]> outList = new ArrayList<>();

        // Trivial cases:
        outList.add(
                new Object[] {
                        new SimpleKeyXsvFuncotationFactory(
                                defaultName,
                                IOUtils.getPath(FuncotatorTestConstants.XSV_CSV_FILE_PATH),
                                "VERSION",
                                ",",
                                0,
                                SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME
                        ),
                        Collections.emptyList(),
                        Collections.emptyList()
                }
        );

        final SimpleKeyXsvFuncotationFactory geneNameXsvFuncotationFactory = new SimpleKeyXsvFuncotationFactory(
                defaultName,
                IOUtils.getPath(FuncotatorTestConstants.XSV_CSV_FILE_PATH),
                "VERSION",
                ",",
                0,
                SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME
        );
        final List<String> emptyFields1 = new ArrayList<>();
        for ( int i = 0 ; i < geneNameXsvFuncotationFactory.getSupportedFuncotationFields().size() ; ++i ) { emptyFields1.add(""); }
        outList.add(
                new Object[] {
                        geneNameXsvFuncotationFactory,
                        Collections.singletonList(
                                new GencodeFuncotationBuilder().setHugoSymbol("NOT THE RIGHT GENE NAME").build()
                        ),
                        Collections.singletonList(
                                TableFuncotation.create(new ArrayList<>(geneNameXsvFuncotationFactory.getSupportedFuncotationFields()), emptyFields1, defaultAltAllele, defaultName, null)
                        )
                }
        );

        final SimpleKeyXsvFuncotationFactory transcriptIdXsvFuncotationFactory = new SimpleKeyXsvFuncotationFactory(
                defaultName,
                IOUtils.getPath(FuncotatorTestConstants.XSV_CSV_FILE_PATH),
                "VERSION",
                ",",
                0,
                SimpleKeyXsvFuncotationFactory.XsvDataKeyType.TRANSCRIPT_ID
        );
        final List<String> emptyFields2 = new ArrayList<>();
        for ( int i = 0 ; i < transcriptIdXsvFuncotationFactory.getSupportedFuncotationFields().size() ; ++i ) { emptyFields2.add(""); }
        outList.add(
            new Object[] {
                    transcriptIdXsvFuncotationFactory,
                    Collections.singletonList(
                            new GencodeFuncotationBuilder().setAnnotationTranscript("NOT THE RIGHT TRANSCRIPT ID").build()
                    ),
                    Collections.singletonList(
                            TableFuncotation.create(new ArrayList<>(transcriptIdXsvFuncotationFactory.getSupportedFuncotationFields()), emptyFields2, defaultAltAllele, defaultName, null)
                    )
            }
        );
        
        // Add in cases from helper function:
        helpPopulateDataForCreateFuncotations(outList, FuncotatorTestConstants.XSV_CSV_FILE_PATH, ",");
        helpPopulateDataForCreateFuncotations(outList, FuncotatorTestConstants.XSV_TSV_FILE_PATH, "\t");
        helpPopulateDataForCreateFuncotations(outList, FuncotatorTestConstants.XSV_DEADBEEFSV_FILE_PATH, "DEADBEEF");
        helpPopulateDataForCreateFuncotations(outList, FuncotatorTestConstants.XSV_PIPESV_FILE_PATH, "|");

        return outList.iterator();
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestGetSupportedFuncotationFields")
    public void testGetSupportedFuncotationFields(final Path inputPath,
                                                  final int headerLinesToIgnore,
                                                  final String name,
                                                  final String delimiter,
                                                  final int keyColumn,
                                                  final SimpleKeyXsvFuncotationFactory.XsvDataKeyType dataKeyType,
                                                  final LinkedHashSet<String> expected) {

        final SimpleKeyXsvFuncotationFactory xsvFuncotationFactory;
        if ( headerLinesToIgnore == 0 ) {
            xsvFuncotationFactory = new SimpleKeyXsvFuncotationFactory(name, inputPath, "VERSION", delimiter, keyColumn, dataKeyType);
        }
        else {
            xsvFuncotationFactory = new SimpleKeyXsvFuncotationFactory(name, inputPath, "VERSION", delimiter, keyColumn, dataKeyType, new LinkedHashMap<>(), headerLinesToIgnore);
        }

        final LinkedHashSet<String> supportedFields = xsvFuncotationFactory.getSupportedFuncotationFields();

        Assert.assertEquals( supportedFields, expected );
    }

    @Test(dataProvider = "provideForGetName")
    public void testGetName( final SimpleKeyXsvFuncotationFactory factory, final String expected ) {
        Assert.assertEquals( factory.getName(), expected );
    }

    @Test(expectedExceptions = GATKException.class)
    public void testCreateFuncotationsNoGencodeInput() {
        final SimpleKeyXsvFuncotationFactory simpleKeyXsvFuncotationFactory =
                new SimpleKeyXsvFuncotationFactory(
                    defaultName,
                    IOUtils.getPath(FuncotatorTestConstants.XSV_CSV_FILE_PATH),
                    "VERSION",
                    ",",
                    0,
                    SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME
                );

        simpleKeyXsvFuncotationFactory.createFuncotationsOnVariant(defaultVariantContext, defaultReferenceContext, Collections.emptyList());
    }

    @Test(dataProvider = "provideDataForCreateFuncotations")
    public void testCreateFuncotations(final SimpleKeyXsvFuncotationFactory xsvFuncotationFactory,
                                       final List<GencodeFuncotation> gencodeFuncotations,
                                       final List<Funcotation> expectedFuncotationsList) {

        final List<Funcotation> funcotations = xsvFuncotationFactory.createFuncotationsOnVariant(defaultVariantContext, defaultReferenceContext, Collections.emptyList(), gencodeFuncotations);

        Assert.assertEquals( funcotations.size(), expectedFuncotationsList.size(),
                "Wrong number of funcotations created (" + funcotations.size() + ")  Expected: " + expectedFuncotationsList.size() );

        for (int i = 0; i < funcotations.size(); ++i) {

            final TableFuncotation computed = (TableFuncotation)funcotations.get(i);
            final TableFuncotation expected = (TableFuncotation)expectedFuncotationsList.get(i);

            Assert.assertEquals(computed.size(), expected.size(),
                    "Funcotations at index " + i + " do not have the same number of elements: computed " + computed.size() + " != " + expected.size() + " expected");

            Assert.assertEquals(computed, expected,
                    "Funcotations at index " + i + " are not equal:" + "\n\tComputed: " + computed + "\n\tExpected: " + expected);

        }
    }

}
