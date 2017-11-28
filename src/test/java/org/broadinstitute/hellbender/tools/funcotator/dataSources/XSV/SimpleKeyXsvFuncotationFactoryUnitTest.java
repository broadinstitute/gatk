package org.broadinstitute.hellbender.tools.funcotator.dataSources.XSV;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
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
    private static final List<List<String>> headerOffsetStrings;

    static {
        headerOffsetStrings = new ArrayList<>(squareSize);
        for ( int i = 0; i < squareSize; ++i ) {
            final List<String> singleRow = new ArrayList<>(squareSize);
            for ( int j = 0; j < squareSize; ++j ) {
                singleRow.add( defaultName + "_R" + (i+1) + "C"+ (j+1) );
            }
            headerOffsetStrings.add(singleRow);
        }
    }

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    private void helpPopulateForGetSupportedFuncotationFields(final List<Object[]> outList,
                                                              final List<Integer> keyColumns,
                                                              final String path,
                                                              final String delim) {
        for ( int i = 0; i < squareSize; ++i ) {
            for ( final int keyCol : keyColumns ) {

                final List<String> headers = new ArrayList<>(headerOffsetStrings.get(i));
                headers.remove(keyCol);

                outList.add(
                        new Object[]{
                                new File(path).toPath(), i, defaultName, delim, keyCol, SimpleKeyXsvFuncotationFactory.XsvDataKeyType.GENE_NAME,
                                new LinkedHashSet<>(headers)
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

        helpPopulateForGetSupportedFuncotationFields(outList, Arrays.asList(1,5,17), FuncotatorTestConstants.XSV_CSV_FILE_PATH, ",");
        helpPopulateForGetSupportedFuncotationFields(outList, Arrays.asList(1,5,17), FuncotatorTestConstants.XSV_TSV_FILE_PATH, "\t");
        helpPopulateForGetSupportedFuncotationFields(outList, Arrays.asList(1,5,17), FuncotatorTestConstants.XSV_DEADBEEFSV_FILE_PATH, "DEADBEEF");

        return outList.iterator();
    }

    //==================================================================================================================
    // Tests:

    // CreateFuncotations (x2)

    // GetSupportedfuncotationFields
    @Test(dataProvider = "provideForTestGetSupportedFuncotationFields")
    public void testGetSupportedfuncotationFields(final Path inputPath,
                                                  final int headerLinesToIgnore,
                                                  final String name,
                                                  final String delimiter,
                                                  final int keyColumn,
                                                  final SimpleKeyXsvFuncotationFactory.XsvDataKeyType dataKeyType,
                                                  final LinkedHashSet<String> expected) {

        final SimpleKeyXsvFuncotationFactory xsvFuncotationFactory =
                new SimpleKeyXsvFuncotationFactory(name, inputPath, delimiter, keyColumn, dataKeyType, headerLinesToIgnore);

        final LinkedHashSet<String> supportedFields = xsvFuncotationFactory.getSupportedFuncotationFields();

        Assert.assertEquals( supportedFields, expected );
    }

    // constructors

}
