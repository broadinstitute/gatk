package org.broadinstitute.hellbender.utils.io;

import com.google.common.collect.ImmutableMap;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.stream.Collectors;

public class ResourceUnitTest extends GATKBaseTest {
    private static final String TEST_SUB_DIR = "org/broadinstitute/hellbender/utils/io/";
    private static final String TEST_RESOURCE_CONFIG = TEST_SUB_DIR + "resource.properties";

    @Test(dataProvider = "propGT")
    public void testGetResourceContentsAsFile(String filePath, Map<String, String> gt) throws IOException {
        // IMPORTANT:  We have seen pathologies where the test passes, but fails when in the jar.
        final File rsc = Resource.getResourceContentsAsFile(filePath);

        final Properties rscProps = new Properties();
        try ( final InputStream inputStream = Files.newInputStream(rsc.toPath(), StandardOpenOption.READ) ) {
            rscProps.load(inputStream);

            final List<String> keys = gt.keySet().stream().sorted().collect(Collectors.toList());
            Assert.assertEquals(keys.stream().map(k -> rscProps.getProperty(k)).collect(Collectors.toList()),
                    keys.stream().map(k -> gt.get(k)).collect(Collectors.toList()));
        }
    }

    @DataProvider(name="propGT")
    public Object[][] createPropGT() {
        return new Object[][]{
                {       TEST_RESOURCE_CONFIG,
                        ImmutableMap.of("key1", "Foo", "key2", "Bar")
                }
        };
    }
}
