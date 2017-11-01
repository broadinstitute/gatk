package org.broadinstitute.hellbender.utils.codecs.sampileup;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;

import static org.broadinstitute.hellbender.utils.BaseUtils.Base.A;
import static org.broadinstitute.hellbender.utils.BaseUtils.Base.C;
import static org.broadinstitute.hellbender.utils.BaseUtils.Base.N;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class SAMPileupFeatureUnitTests extends GATKBaseTest {

    @Test
    public void testGetPileupString() {
        final SAMPileupElement element1 = new SAMPileupElement(A.base, (byte) 40);
        final SAMPileupElement element2 = new SAMPileupElement(C.base, (byte) 26);
        Assert.assertEquals(new SAMPileupFeature("chr1", 100, C.base, Collections.singletonList(element1)).getPileupString(), "chr1 100 C A I");
        Assert.assertEquals(new SAMPileupFeature("chr2", 200, A.base, Collections.singletonList(element2)).getPileupString(), "chr2 200 A C ;");
        Assert.assertEquals(new SAMPileupFeature("chr3", 1000, N.base, Arrays.asList(element1, element2)).getPileupString(), "chr3 1000 N AC I;");
    }

}
