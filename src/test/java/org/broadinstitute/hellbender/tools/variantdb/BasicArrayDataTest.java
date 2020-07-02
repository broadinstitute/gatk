package org.broadinstitute.hellbender.tools.variantdb;

import org.testng.annotations.Test;
import static org.testng.Assert.*;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.variantdb.BasicArrayData.ArrayGenotype;

public final class BasicArrayDataTest {

    @Test
    public void testEncodeDecode() {
        BasicArrayData original = new BasicArrayData(12345, 67890, ArrayGenotype.AB);
        long bits = original.encode();        
        BasicArrayData r = new BasicArrayData(bits);

        assertEquals(r.genotype, original.genotype);
        assertEquals(r.probeId, original.probeId);
        assertEquals(r.sampleId, original.sampleId);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testSampleIdTooSmall() {
        BasicArrayData original = new BasicArrayData(-1, 67890, ArrayGenotype.AB);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testSampleIdTooBig() {
        BasicArrayData original = new BasicArrayData(BasicArrayData.MAX_SAMPLE_ID_VALUE+1, 67890, ArrayGenotype.AB);
    }

    @Test
    public void testSampleIdAtLimits() {
        BasicArrayData original = new BasicArrayData(0, 67890, ArrayGenotype.AB);
        BasicArrayData r = new BasicArrayData(original.encode());
        assertEquals(r.sampleId, original.sampleId);

        original = new BasicArrayData(BasicArrayData.MAX_SAMPLE_ID_VALUE, 67890, ArrayGenotype.AB);
        r = new BasicArrayData(original.encode());
        assertEquals(r.sampleId, original.sampleId);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testProbeIdTooSmall() {
        BasicArrayData original = new BasicArrayData(1, -1, ArrayGenotype.AB);
    }

    @Test(expectedExceptions = GATKException.class)
    public void testProbeIdTooBig() {
        BasicArrayData original = new BasicArrayData(1, BasicArrayData.MAX_PROBE_ID_VALUE+1, ArrayGenotype.AB);
    }

    @Test
    public void testProbeIdAtLimits() {
        BasicArrayData original = new BasicArrayData(1, 0, ArrayGenotype.AB);
        BasicArrayData r = new BasicArrayData(original.encode());
        assertEquals(r.probeId, original.probeId);

        original = new BasicArrayData(1, BasicArrayData.MAX_PROBE_ID_VALUE, ArrayGenotype.AB);
        r = new BasicArrayData(original.encode());
        assertEquals(r.probeId, original.probeId);
    }
}
