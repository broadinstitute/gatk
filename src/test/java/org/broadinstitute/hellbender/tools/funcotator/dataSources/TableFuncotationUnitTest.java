package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A Unit Test class for {@link TableFuncotation}
 * Created by jonn on 11/28/17.
 */
public class TableFuncotationUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideForTestGetAltAllele() {
        return new Object[][] {
                {
                        TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("A", false), "TableFuncotationUnitTest", null),
                        Allele.create("A", false),
                        true
                },
                {
                        TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("C", false), "TableFuncotationUnitTest", null),
                        Allele.create("C", false),
                        true
                },
                {
                        TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("G", false), "TableFuncotationUnitTest", null),
                        Allele.create("G", false),
                        true
                },
                {
                        TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("T", false), "TableFuncotationUnitTest", null),
                        Allele.create("T", false),
                        true
                },
                {
                        TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("C", false), "TableFuncotationUnitTest", null),
                        Allele.create("A", false),
                        false
                },
                {
                        TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("G", false), "TableFuncotationUnitTest", null),
                        Allele.create("C", false),
                        false
                },
                {
                        TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("T", false), "TableFuncotationUnitTest", null),
                        Allele.create("G", false),
                        false
                },
                {
                        TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("G", true), "TableFuncotationUnitTest", null),
                        Allele.create("G", false),
                        false
                },
                {
                        TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("T", true), "TableFuncotationUnitTest", null),
                        Allele.create("T", false),
                        false
                },
        };
    }

    @DataProvider
    Object[][] provideForTestGet() {
        return new Object[][]{
                { "X" },
                { "A" },
                { "B" },
                { "Y" },
                { "L" },
                { "R" },
                { "Select" },
                { "Start" },
        };
    }

    @DataProvider
    Object[][] provideDataForTestSetFieldSerializationOverrideValue() {

        final TableFuncotation funcotation =
                TableFuncotation.create(
                    Arrays.asList("A", "B", "C"), Arrays.asList("1", "2", "3"), Allele.create("A", false), "TableFuncotation", null
                );

        return new Object[][] {
                { funcotation, "A", null },
                { funcotation, "A", "A" },
                { funcotation, "A", "B" },
                { funcotation, "A", "C" },
                { funcotation, "B", "A" },
                { funcotation, "B", "B" },
                { funcotation, "B", "C" },
                { funcotation, "C", "A" },
                { funcotation, "C", "B" },
                { funcotation, "C", "C" },
        };
    }

    @DataProvider
    Object[][] provideForTestSerializeToVcfString() {
        return new Object[][] {
                {
                    TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("A", false), "Empty", null),
                    ""
                },
                {
                    TableFuncotation.create(Collections.singletonList("A"), Collections.singletonList(("1")), Allele.create("A", false), "OneVal", null),
                    "1"
                },
                {
                    TableFuncotation.create(Arrays.asList("A", "B"), Arrays.asList("1", "2"), Allele.create("A", false), "TwoVals", null),
                    "1|2"
                },
                {
                    TableFuncotation.create(Arrays.asList("A", "B", "C"), Arrays.asList("1", "2", "3"), Allele.create("A", false), "ThreeVals", null),
                    "1|2|3"
                },
        };
    }

    @DataProvider
    Object[][] provideListOfStrings() {
        return new Object[][] {
                { Collections.singletonList("A") },
                { Arrays.asList("A", "B") },
                { Arrays.asList("A", "B", "C") },
                { Arrays.asList("A", "B", "C", "D") },
                { Arrays.asList("A", "B", "C", "D", "E") },
        };
    }

    @DataProvider
    Object[][] provideForTestGetFieldNames() {
        //final TableFuncotation tableFuncotation, final LinkedHashSet<String> expected
        return new Object[][] {
                {
                    TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("A", false), "Empty", null),
                    new LinkedHashSet<>(Collections.emptyList())
                },
                {
                    TableFuncotation.create(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), Allele.create("A", false), "OneField", null),
                    new LinkedHashSet<>(Collections.singletonList("TESTFIELD"))
                },
                {
                    TableFuncotation.create(Arrays.asList("TESTFIELD1", "TESTFIELD2"), Arrays.asList("TESTVAL1", "TESTVAL2"), Allele.create("A", false), "TwoFields", null),
                    new LinkedHashSet<>(Arrays.asList("TESTFIELD1", "TESTFIELD2"))
                },
        };
    }

    @DataProvider
    Object[][] provideForTestGetField() {
        //final TableFuncotation tableFuncotation, final String fieldName, final String expected
        return new Object[][] {
                {
                    TableFuncotation.create(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), Allele.create("A", false), "OneField", null),
                    "TESTFIELD",
                    "TESTVAL"
                },
                {
                    TableFuncotation.create(Arrays.asList("TESTFIELD1", "TESTFIELD2"), Arrays.asList("TESTVAL1", "TESTVAL2"), Allele.create("A", false), "TwoFields", null),
                    "TESTFIELD1",
                    "TESTVAL1"
                },
                {
                    TableFuncotation.create(Arrays.asList("TESTFIELD1", "TESTFIELD2"), Arrays.asList("TESTVAL1", "TESTVAL2"), Allele.create("A", false), "TwoFields", null),
                    "TESTFIELD2",
                    "TESTVAL2"
                },
                {
                    TableFuncotation.create(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), Allele.create("A", false), "ThreeFields", null),
                    "TESTFIELD1",
                    "TESTVAL1"
                },
                {
                    TableFuncotation.create(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), Allele.create("A", false), "ThreeFields", null),
                    "TESTFIELD2",
                    "TESTVAL2"
                },
                {
                    TableFuncotation.create(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), Allele.create("A", false), "ThreeFields", null),
                    "TESTFIELD3",
                    "TESTVAL3"
                },
        };
    }

    @DataProvider
    Object[][] provideForTestGetFieldFail() {
        //final TableFuncotation tableFuncotation, final String fieldName, final String expected
        return new Object[][] {
                {
                    TableFuncotation.create(Collections.emptyList(), Collections.emptyList(), Allele.create("A", false), "Empty", null),
                    "TESTFIELD_OMICRON"
                },
                {
                    TableFuncotation.create(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), Allele.create("A", false), "OneField", null),
                    "testfield"
                },
                {
                    TableFuncotation.create(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), Allele.create("A", false), "OneField", null),
                    "table_TESTFIELD"
                },
                {
                    TableFuncotation.create(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), Allele.create("A", false), "ThreeFields", null),
                    "TESTFIELD4"
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestGetAltAllele")
    public void testGetAltAllele( final TableFuncotation funcotation, final Allele altAllele, final boolean expected) {
        Assert.assertEquals( funcotation.getAltAllele().equals(altAllele), expected );
    }

    @Test(dataProvider = "provideForTestGet")
    public void testGet(final String fieldValue) {

        final String fieldName = "PLACEHOLDER";

        final TableFuncotation funcotation = TableFuncotation.create( Collections.singletonList(fieldName), Collections.singletonList(fieldValue), Allele.create("A", false), fieldName , null);
        Assert.assertEquals( funcotation.get(fieldName), fieldValue );
    }

    @Test(dataProvider = "provideDataForTestSetFieldSerializationOverrideValue")
    public void testSetFieldSerializationOverrideValue(final TableFuncotation funcotation,
                                                       final String fieldName,
                                                       String overrideValue) {
        if ( overrideValue == null ) {
            overrideValue = funcotation.get(fieldName);
        }
        else {
            funcotation.setFieldSerializationOverrideValue(fieldName, overrideValue);
        }
        Assert.assertEquals( funcotation.get(fieldName), overrideValue );
    }

    @Test(dataProvider = "provideForTestSerializeToVcfString")
    public void testSerializeToVcfString(final TableFuncotation funcotation,
                                         final String expected) {
        Assert.assertEquals( funcotation.serializeToVcfString(), expected );
    }

    @Test(dataProvider = "provideListOfStrings")
    public void testKeySet(final List<String> kvNames) {
        final TableFuncotation funcotation = TableFuncotation.create(kvNames, kvNames.stream().map(s -> s + "VVV").collect(Collectors.toList()), Allele.create("A", false), "ListFuncotation", null);
        Assert.assertEquals( funcotation.keySet(), kvNames);
    }

    @Test(dataProvider = "provideListOfStrings")
    public void testValues(final List<String> kvNames) {
        final TableFuncotation funcotation = TableFuncotation.create(kvNames.stream().map(s -> s + "KEY").collect(Collectors.toList()), kvNames, Allele.create("A", false), "ListFuncotation", null);
        Assert.assertEquals( funcotation.values(), kvNames);
    }

    @Test(dataProvider = "provideForTestGetFieldNames")
    public void testGetFieldNames(final TableFuncotation tableFuncotation, final LinkedHashSet<String> expected) {
        Assert.assertEquals(tableFuncotation.getFieldNames(), expected);
    }

    @Test(dataProvider = "provideForTestGetField")
    public void testGetField(final TableFuncotation tableFuncotation, final String fieldName, final String expected) {
        Assert.assertEquals(tableFuncotation.getField(fieldName), expected);
    }

    @Test(dataProvider = "provideForTestGetFieldFail", expectedExceptions = GATKException.class)
    public void testGetFieldFail(final TableFuncotation tableFuncotation, final String fieldName) {
        tableFuncotation.getField(fieldName);
    }

}
