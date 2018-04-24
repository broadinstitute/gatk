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
                        new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("A", false), "TableFuncotationUnitTest"),
                        Allele.create("A", false),
                        true
                },
                {
                        new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("C", false), "TableFuncotationUnitTest"),
                        Allele.create("C", false),
                        true
                },
                {
                        new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("G", false), "TableFuncotationUnitTest"),
                        Allele.create("G", false),
                        true
                },
                {
                        new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("T", false), "TableFuncotationUnitTest"),
                        Allele.create("T", false),
                        true
                },
                {
                        new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("C", false), "TableFuncotationUnitTest"),
                        Allele.create("A", false),
                        false
                },
                {
                        new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("G", false), "TableFuncotationUnitTest"),
                        Allele.create("C", false),
                        false
                },
                {
                        new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("T", false), "TableFuncotationUnitTest"),
                        Allele.create("G", false),
                        false
                },
                {
                        new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("G", true), "TableFuncotationUnitTest"),
                        Allele.create("G", false),
                        false
                },
                {
                        new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("T", true), "TableFuncotationUnitTest"),
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
                new TableFuncotation(
                    Arrays.asList("A", "B", "C"), Arrays.asList("1", "2", "3"), Allele.create("A", false), "TableFuncotation"
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
                    new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("A", false), "Empty"),
                    ""
                },
                {
                    new TableFuncotation(Collections.singletonList("A"), Collections.singletonList(("1")), Allele.create("A", false), "OneVal"),
                    "1"
                },
                {
                    new TableFuncotation(Arrays.asList("A", "B"), Arrays.asList("1", "2"), Allele.create("A", false), "TwoVals"),
                    "1|2"
                },
                {
                    new TableFuncotation(Arrays.asList("A", "B", "C"), Arrays.asList("1", "2", "3"), Allele.create("A", false), "ThreeVals"),
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
                    new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("A", false), "Empty"),
                    new LinkedHashSet<>(Collections.emptyList())
                },
                {
                    new TableFuncotation(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), Allele.create("A", false), "OneField"),
                    new LinkedHashSet<>(Collections.singletonList("TESTFIELD"))
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2"), Arrays.asList("TESTVAL1", "TESTVAL2"), Allele.create("A", false), "TwoFields"),
                    new LinkedHashSet<>(Arrays.asList("TESTFIELD1", "TESTFIELD2"))
                },
        };
    }

    @DataProvider
    Object[][] provideForTestGetField() {
        //final TableFuncotation tableFuncotation, final String fieldName, final String expected
        return new Object[][] {
                {
                    new TableFuncotation(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), Allele.create("A", false), "OneField"),
                    "TESTFIELD",
                    "TESTVAL"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2"), Arrays.asList("TESTVAL1", "TESTVAL2"), Allele.create("A", false), "TwoFields"),
                    "TESTFIELD1",
                    "TESTVAL1"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2"), Arrays.asList("TESTVAL1", "TESTVAL2"), Allele.create("A", false), "TwoFields"),
                    "TESTFIELD2",
                    "TESTVAL2"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), Allele.create("A", false), "ThreeFields"),
                    "TESTFIELD1",
                    "TESTVAL1"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), Allele.create("A", false), "ThreeFields"),
                    "TESTFIELD2",
                    "TESTVAL2"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), Allele.create("A", false), "ThreeFields"),
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
                    new TableFuncotation(Collections.emptyList(), Collections.emptyList(), Allele.create("A", false), "Empty"),
                    "TESTFIELD_OMICRON"
                },
                {
                    new TableFuncotation(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), Allele.create("A", false), "OneField"),
                    "testfield"
                },
                {
                    new TableFuncotation(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), Allele.create("A", false), "OneField"),
                    "table_TESTFIELD"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), Allele.create("A", false), "ThreeFields"),
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

        final TableFuncotation funcotation = new TableFuncotation( Collections.singletonList(fieldName), Collections.singletonList(fieldValue), Allele.create("A", false), fieldName );
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
        final TableFuncotation funcotation = new TableFuncotation(kvNames, kvNames.stream().map(s -> s + "VVV").collect(Collectors.toList()), Allele.create("A", false), "ListFuncotation");
        Assert.assertEquals( funcotation.keySet(), kvNames);
    }

    @Test(dataProvider = "provideListOfStrings")
    public void testValues(final List<String> kvNames) {
        final TableFuncotation funcotation = new TableFuncotation(kvNames.stream().map(s -> s + "KEY").collect(Collectors.toList()), kvNames, Allele.create("A", false), "ListFuncotation");
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
