package org.broadinstitute.hellbender.engine.transformers;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.LinkedHashMap;
import java.util.Map;

public class IdentityTransformerUnitTest extends GATKBaseTest {

    final static String ATTR1 = "A1";
    final static String ATTR2 = "A2";

    @Test
    public void testIdentity() {

        // create a read w/ arbitrary attributes
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "Zuul", 0, 2,2);
        read.setAttribute("CN", getClass().getName());
        read.setAttribute("PN", getClass().getPackage().getName());

        // snapshot bases, quality and attributes (no way to get all attributes on GATKRead!)
        final String bases = read.getBasesString();
        final String quals = new String(read.getBaseQualitiesNoCopy());
        final String a1 = read.getAttributeAsString(ATTR1);
        final String a2 = read.getAttributeAsString(ATTR2);

        // transform
        final ReadTransformer readTransformer = new IdentityTransformer();
        final GATKRead dst = readTransformer.apply(read);

        // check
        Assert.assertEquals(dst.getBasesString(), bases);
        Assert.assertEquals(new String(dst.getBaseQualitiesNoCopy()), quals);
        Assert.assertEquals(dst.getAttributeAsString(ATTR1), a1);
        Assert.assertEquals(dst.getAttributeAsString(ATTR2), a2);
    }

}
