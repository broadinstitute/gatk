package org.broadinstitute.hellbender.tools.walkers.misc;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramUnitTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.junit.Assert;
import org.testng.annotations.Test;

import java.nio.file.Paths;

/**
 * Created by tsato on 1/25/18.
 */
public class CreateHomopolymerIntervalsUnitTest extends CommandLineProgramUnitTest{

    @Test
    public void testGetTrailingBases(){
        final int numberOfBases = 5;
        final SimpleInterval testInterval = new SimpleInterval("20", 100_000, 100_000);
        final ReferenceContext referenceContext = new ReferenceContext(new ReferenceFileSource(Paths.get(b37_reference_20_21)), testInterval);
        final byte[] trailingBases = CreateHomopolymerIntervals.getTrailingBases(referenceContext, numberOfBases);
        Assert.assertArrayEquals(trailingBases, new byte[]{ 'C', 'T', 'T', 'T', 'C' });

    }

    @Test
    public void testGetPrecedingBases(){
        final int numberOfBases = 5;
        final SimpleInterval testInterval = new SimpleInterval("20", 100_000, 100_000);
        final ReferenceContext referenceContext = new ReferenceContext(new ReferenceFileSource(Paths.get(b37_reference_20_21)), testInterval);
        final byte[] precedingBases = CreateHomopolymerIntervals.getPrecedingBases(referenceContext, numberOfBases);
        Assert.assertArrayEquals(precedingBases, new byte[]{ 'A', 'T', 'T', 'T', 'T' });

    }

}