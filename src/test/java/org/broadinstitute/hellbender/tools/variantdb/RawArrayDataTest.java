package org.broadinstitute.hellbender.tools.variantdb;

//import org.testng.annotations.Test;
import org.testng.annotations.Test;
import static org.testng.Assert.*;

import org.broadinstitute.hellbender.tools.variantdb.RawArrayData.ArrayGenotype;

public final class RawArrayDataTest {

    @Test
    public void testEncodeDecode() {
        RawArrayData original = new RawArrayData();
        float e = 1.0f / 256.0f;

        original.lrr = e;
        original.baf = 2*e;
        original.normx = 3*e;
        original.normy = 4*e;
        original.genotype = ArrayGenotype.AB;
        original.probeId = 123456;

        long bits = original.encode();        
        RawArrayData r = RawArrayData.decode(bits);

        assertEquals(r.lrr, original.lrr);
        assertEquals(r.baf, original.baf);
        assertEquals(r.normx, original.normx);
        assertEquals(r.normy, original.normy);        
        assertEquals(r.genotype, original.genotype);
        assertEquals(r.probeId, original.probeId);
    }
}
