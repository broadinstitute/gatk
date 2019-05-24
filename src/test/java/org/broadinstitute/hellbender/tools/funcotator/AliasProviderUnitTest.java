package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadataUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;

public class AliasProviderUnitTest extends GATKBaseTest {

    @DataProvider
    public Object[][] provideColumnNameToFieldName() {
        return new Object[][] {

                {
                    // Everything is empty
                    new LinkedHashMap<>(), FuncotationMap.createNoTranscriptInfo(Collections.emptyList()), FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY,
                        new LinkedHashMap<>()
                },{
                    // Trivial alias works
                    FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("FOO"), Collections.singletonList(Arrays.asList("BAR", "BAZ"))),
                    FuncotationMap.createNoTranscriptInfo(
                            Collections.singletonList(
                                    createDummyTestFuncotations().get(0)
                            )
                    ),
                    FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY,
                    FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("FOO"), Collections.singletonList("BAZ"))
                },{
                    // Trivial alias works with multiple funcotations
                    FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("FOO"), Collections.singletonList(Arrays.asList("BAR", "BAZ", "BAZ2"))),
                    FuncotationMap.createNoTranscriptInfo(
                            createDummyTestFuncotations()
                    ),
                    FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY,
                    FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("FOO"), Collections.singletonList("BAZ"))
                },{
                    // Trivial alias works with multiple funcotations
                    FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("FOO"), Collections.singletonList(Arrays.asList("BAR", "BAZ2", "BAZ2"))),
                    FuncotationMap.createNoTranscriptInfo(
                            createDummyTestFuncotations()
                    ),
                    FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY,
                    FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("FOO"), Collections.singletonList("BAZ2"))
                },{
                    // No alias found (returned alias should be empty string)
                    FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("FOO"), Collections.singletonList(Arrays.asList("XXX", "YYY", "ZZZ"))),
                    FuncotationMap.createNoTranscriptInfo(
                            createDummyTestFuncotations()
                    ),
                    FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY,
                    FuncotatorUtils.createLinkedHashMapFromLists(Collections.singletonList("FOO"), Collections.singletonList(""))
                },{
                    // Multiple aliases
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("FIELD1", "FIELD2"), Arrays.asList(Arrays.asList("XXX", "GUZ2", "ZZZ"), Arrays.asList("BAZ", "BBB", "CCC"))),
                    FuncotationMap.createNoTranscriptInfo(
                            createDummyTestFuncotations()
                    ),
                    FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY,
                    FuncotatorUtils.createLinkedHashMapFromLists(Arrays.asList("FIELD1", "FIELD2"), Arrays.asList("GUZ2", "BAZ"))
                }
        };
    }

    private List<Funcotation> createDummyTestFuncotations() {
        return Arrays.asList(
            TableFuncotation.create(Arrays.asList("BAZ", "GUZ"), Arrays.asList("VAL_BAZ", "VAL_GUZ"), Allele.create("T"), "TEST",
                FuncotationMetadataUtils.createWithUnknownAttributes(Arrays.asList("BAZ", "GUZ"))),
            TableFuncotation.create(Arrays.asList("BAZ2", "GUZ2"), Arrays.asList("VAL_BAZ2", "VAL_GUZ2"), Allele.create("T"), "TEST2",
                FuncotationMetadataUtils.createWithUnknownAttributes(Arrays.asList("BAZ2", "GUZ2")))
        );
    }

    @Test(dataProvider = "provideColumnNameToFieldName")
    public void testColumnNameToFieldName(final LinkedHashMap<String, List<String>> aliasMap, final FuncotationMap funcotationMap, final String txId, final LinkedHashMap<String, String> gt) {
        final AliasProvider aliasProvider = new AliasProvider(aliasMap);
        final LinkedHashMap<String, String> guess =
                aliasProvider.createColumnNameToFieldNameMap(funcotationMap, txId);
        Assert.assertEquals(guess, gt);
        Assert.assertEquals(aliasProvider.getFields(), aliasMap.keySet());
    }
}
