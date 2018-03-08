package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.collections.MapUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.tools.funcotator.TranscriptSelectionMode;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.*;

/**
 * Unit test class for the {@link MafOutputRenderer}.
 * Created by jonn on 1/22/18.
 */
public class MafOutputRendererUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    private MafOutputRenderer createMafOutputRenderer() {

        final Map<Path, Properties> configData =
                DataSourceUtils.getAndValidateDataSourcesFromPaths(FuncotatorTestConstants.REFERENCE_VERSION_HG19, Collections.singletonList(FuncotatorTestConstants.FUNCOTATOR_DATA_SOURCES_MAIN_FOLDER));

        return new MafOutputRenderer(
                IOUtils.getPath(getSafeNonExistentFile("TestMafOutputFile").toURI().toString()),
                DataSourceUtils.createDataSourceFuncotationFactoriesForDataSources(configData, new LinkedHashMap<>(), TranscriptSelectionMode.BEST_EFFECT, new HashSet<>()),
                new VCFHeader(),
                new LinkedHashMap<>(),
                new LinkedHashMap<>(),
                new HashSet<>()
        );
    }

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Object[][] provideForTestReplaceFuncotationValuesWithMafCompliantValues() {
        return new Object[][] {
                // Empty maps:
                {
                    new HashMap<>(),
                    new LinkedHashMap<>()
                },
                // Singleton map that doesn't contain a replaceable element:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        {"", ""},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"", ""},
                                }
                        )
                },
                // Singleton map containing a replaced value:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        {"", ""},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"", ""},
                                }
                        )
                },
                // Map with multiple elements, none of which are replaced:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        {"", ""},
                                        {"", ""},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"", ""},
                                }
                        )
                },
                // Map with multiple elements, some of which are replaced:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        {"", ""},
                                        {"", ""},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"", ""},
                                }
                        )
                },
                // Map with multiple elements, all of which are replaced:
                {
                        MapUtils.putAll(new HashMap<String, Object>(),
                                new Object[][] {
                                        {"", ""},
                                        {"", ""},
                                }
                        ),
                        MapUtils.putAll(new LinkedHashMap<String, String>(),
                                new Object[][] {
                                        {"", ""},
                                }
                        )
                }
        };
    }

    @DataProvider
    private Object[][] provideForMafTransform() {
        return new Object[][] {
                // Empty Strings:
                {
                    "", "", ""
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestReplaceFuncotationValuesWithMafCompliantValues")
    public void testReplaceFuncotationValuesWithMafCompliantValues(final Map<String, Object> outputMap, final LinkedHashMap<String, String> expected ) {
        final LinkedHashMap<String, String> compliantMap = createMafOutputRenderer().replaceFuncotationValuesWithMafCompliantValues(outputMap);
        Assert.assertEquals(compliantMap, expected);
    }

    @Test(dataProvider = "provideForMafTransform")
    public void testMafTransform(final String key, final String value, final String expectedValue) {
        final String transformedValue = createMafOutputRenderer().mafTransform(key, value);
        Assert.assertEquals(transformedValue, expectedValue);
    }
}
