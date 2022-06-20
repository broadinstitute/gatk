package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class ByteMapperUnitTest extends GATKBaseTest {

    final public String PACKAGE_FOLDER = getClass().getPackage().getName().replace('.', '/');
    final public String TEST_FILE = PACKAGE_FOLDER + "/ByteMapper_test_mapping";
    final public String TEST_FILE_NO_DEFAULT = PACKAGE_FOLDER + "/ByteMapper_test_mapping_no_default";

    @DataProvider(name = "mapping")
    public Object[][] testMappingDataProvider() {
        return new Object[][] {
                {
                    new GATKPath(publicTestDir + "/" + TEST_FILE),
                        new Byte[] { 0, 1, 2, 3, 4, 5 },
                        new Byte[] { 0, 0, 0, 3, 3, 10}
                },
               {
                   new GATKPath(publicTestDir + "/" + TEST_FILE_NO_DEFAULT),
                        new Byte[] { 0, 1, 2, 3, 4, 5 },
                        new Byte[] { 0, 0, 0, 3, 3, 5}
                },
        };
    }

    @Test(dataProvider = "mapping")
    public void testMapping(final GATKPath path, final Byte[] src, final Byte[] dst) {

        final ByteMapper      mapper = new ByteMapper(path);
        Utils.validate(src.length == dst.length, "src and dst should be the same length");
        for ( int i = 0 ; i < src.length ; i++ ) {
            final Byte b = mapper.map(src[i]);
            Assert.assertEquals(b, dst[i]);
        }
    }
}
