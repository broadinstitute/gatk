package org.broadinstitute.hellbender.tools.variantdb.arrays;

import org.testng.annotations.Test;
import static org.testng.Assert.*;

public final class RawArrayDataTest {

    @Test
    public void testEncodeDecode() {
        RawArrayData original = new RawArrayData(1.059f, 2.849f, 0.988f, -0.035f);
        long bits = original.encode();
        
        RawArrayData r = new RawArrayData(bits);
        assertEquals(r.normx, original.normx, 0.001);
        assertEquals(r.normy, original.normy, 0.001);
        assertEquals(r.lrr, original.lrr, 0.001);
        assertEquals(r.baf, original.baf, 0.001);
    }

    @Test
    public void testEncodeDecodeNulls() {
        RawArrayData original = new RawArrayData(0.0f, 0.0f, null, null); 
        long bits = original.encode();
        
        System.out.println("BITS:" + bits);
        
        RawArrayData r = new RawArrayData(bits);
        assertEquals(r.normx, original.normx, 0.001);
        assertEquals(r.normy, original.normy, 0.001);
        assertNull(r.lrr);
        assertNull(r.baf);
    }
}
