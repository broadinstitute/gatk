package org.broadinstitute.hellbender.tools.gvs.utils;

import org.broadinstitute.hellbender.tools.gvs.common.ChromosomeEnum;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class SchemaUtilsTest {

    private final String CHROMOSOME = "1";
    private final int POSITION = 1500;
    private final int LENGTH = 500;
    private final int STATE = 5;

    /*
    The expected encoding of the values above
    top 16 bits (chromosome):    00000000 00000001
    next 32 bits (position):     00000000 00000000 00000101 11011100
    next 12 bits (length):       0001 11110100
    last 4 bits (state):         0101
     */
    private final long EXPECTED_COMPRESSED_VALUE_REAL_VALUES = 0b0000000000000001000000000000000000000101110111000001111101000101L;

    /*
    There isn't a chromosome "0," but we'll zero out all the other values and ensure that they pack that way
     */
    private final long EXPECTED_COMPRESSED_VALUE_ZERO_INPUTS = 0b0000000000000001000000000000000000000000000000000000000000000000L;

    @BeforeClass
    public void setup() {
        // Must set the ref version for the SchemaUtils code to function correctly.
        ChromosomeEnum.setRefVersion("37");
    }

    @Test
    public void testCompressedReferences() {
        long compressedRefValue = SchemaUtils.encodeCompressedRefBlock(CHROMOSOME, POSITION, LENGTH, STATE);
        Assert.assertEquals(compressedRefValue, EXPECTED_COMPRESSED_VALUE_REAL_VALUES);
    }

    @Test
    public void testCompressedReferencesZeroInputs() {
        long compressedRefValue = SchemaUtils.encodeCompressedRefBlock("1", 0, 0, 0);
        Assert.assertEquals(compressedRefValue, EXPECTED_COMPRESSED_VALUE_ZERO_INPUTS);
    }
}
