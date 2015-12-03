package org.broadinstitute.hellbender.utils.baq;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class BAQUnitTest extends BaseTest{

    @Test //regression test for https://github.com/broadinstitute/gatk/issues/1234
    public void testGetReferenceWindowForReadStopsAtContigStart() throws Exception {
        final boolean includeClippedBases = false;
        final int bandwidth= 7;
        GATKRead read = ArtificialReadUtils.createArtificialRead("10M");

        final int start = 2;
        final Locatable pos = new SimpleInterval("1", start, start+read.getLength());
        read.setPosition(pos);
        final SimpleInterval referenceWindowForRead = BAQ.getReferenceWindowForRead(read, bandwidth, includeClippedBases);
        SimpleInterval refWindow = new SimpleInterval("1", 1, read.getEnd() + (bandwidth/2));  //start is at 1 because we hit the front end of the reference
        Assert.assertEquals(referenceWindowForRead, refWindow);
    }
}
