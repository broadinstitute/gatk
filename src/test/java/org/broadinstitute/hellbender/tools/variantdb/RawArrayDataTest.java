package org.broadinstitute.hellbender.tools.variantdb;

//import org.testng.annotations.Test;
import org.testng.annotations.Test;
import static org.testng.Assert.*;

import org.broadinstitute.hellbender.tools.variantdb.RawArrayData.ArrayGenotype;
import static org.broadinstitute.hellbender.tools.variantdb.BinaryUtils.*;
import static org.broadinstitute.hellbender.tools.variantdb.RawArrayData.*;


public final class RawArrayDataTest {

    @Test
    public void testEncodeDecode() {
        RawArrayData original = new RawArrayData();

        original.lrr = -0.035f;
        original.baf = 0.988f;
        original.normx = 0.059f;
        original.normy = 0.849f;
        original.genotype = ArrayGenotype.AA;
        original.probeId = 222324;

        long bits = original.encode();        
        RawArrayData r = RawArrayData.decode(bits);

        assertEquals(r.lrr, original.lrr, ( (LRR_MAX - LRR_MIN) / 255.0 ));
        assertEquals(r.baf, original.baf, ( (BAF_MAX - BAF_MIN) / 255.0 ));
        assertEquals(r.normx, original.normx, ( (NORMX_MAX - NORMX_MIN) / 255.0 ));
        assertEquals(r.normy, original.normy, ( (NORMY_MAX - NORMY_MIN) / 255.0 ));
        assertEquals(r.genotype, original.genotype);
        assertEquals(r.probeId, original.probeId);
    }

    @Test
    public void testBasicFloat() { 
        float orig = (1f/255f);
        float f = decodeFrom8Bits((int) encodeTo8Bits(orig, 0, 1), 0, 1);
        assertEquals(orig, f, ( (1 - 0) / 255.0 ));
    }

    @Test
    public void testBasicNull() { 
        Float orig = null;
        Float f = decodeFrom8Bits((int) encodeTo8Bits(orig, 0, 1), 0, 1);
        assertNull(f);
    }

    @Test
    public void testBasicMax() { 
        Float orig = 1.0f;
        Float f = decodeFrom8Bits((int) encodeTo8Bits(orig, 0, 1), 0, 1);
        assertEquals(orig, f);
    }

    @Test
    public void testBasicMin() { 
        Float orig = 0.0f;
        Float f = decodeFrom8Bits((int) encodeTo8Bits(orig, 0, 1), 0, 1);
        assertEquals(orig, f);
    }

    @Test
    public void testBasicLrr() { 
        Float orig = 0.058f;        
        Float f = decodeFrom8Bits((int) encodeTo8Bits(orig, LRR_MIN, LRR_MAX), LRR_MIN, LRR_MAX);
        assertEquals(orig, f, ( (LRR_MAX - LRR_MIN) / 255.0 ));
    }
}
