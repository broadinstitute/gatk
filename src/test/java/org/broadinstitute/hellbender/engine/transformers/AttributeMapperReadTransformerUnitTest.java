package org.broadinstitute.hellbender.engine.transformers;

import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class AttributeMapperReadTransformerUnitTest extends GATKBaseTest {

    final static String ATTR1 = "A1";
    final static String ATTR2 = "A2";

    final public String PACKAGE_FOLDER = getClass().getPackage().getName().replace('.', '/');
    final public String TEST_FILE = PACKAGE_FOLDER + "/test_mapping";

    @DataProvider(name = "mapping")
    public Object[][] testMappingDataProvider() {
        return new Object[][] {
                {
                        new GATKPath(publicTestDir + "/" + TEST_FILE),
                        new Byte[] { 0, 1, 2, 3, 4, 5 },
                        new Byte[] { 0, 0, 0, 3, 3, 10}
                }
        };
    }

    @Test(dataProvider = "mapping")
    public void testQualityMapping(final GATKPath path, final Byte[] src, final Byte[] dst) {


        // create a read w/ attributes
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "Zuul", 0, 2,src.length);
        read.setAttribute(ATTR1, ArrayUtils.toPrimitive(src));
        read.setAttribute(ATTR2, org.bouncycastle.util.Arrays.reverse(ArrayUtils.toPrimitive(src)));

        // transform
        final AttributeMapperReadTransformer readTransformer = new AttributeMapperReadTransformer();
        readTransformer.attributeMappingFile = Arrays.asList(path, path);
        readTransformer.attributeMappingName = Arrays.asList(ATTR1, ATTR2);
        final GATKRead read2 = readTransformer.apply(read);

        // check
        Assert.assertEquals(read2.getAttributeAsByteArray(ATTR1), ArrayUtils.toPrimitive(dst));
        Assert.assertEquals(read2.getAttributeAsByteArray(ATTR2), org.bouncycastle.util.Arrays.reverse(ArrayUtils.toPrimitive(dst)));
    }
}
