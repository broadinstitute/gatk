package org.broadinstitute.hellbender.tools.walkers.fasta;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class CountBasesInReferenceIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testCountBasesInReferenceWithInterval(){
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg19MiniReference))
                .add("L", "1:5000-6000");
        final CountBasesInReference walker = new CountBasesInReference();
        walker.instanceMain(args.getArgsArray());
        final long[] baseCounts = walker.baseCounts;

        Assert.assertEquals(baseCounts['N'], 1001L);
    }

    @Test
    public void testCountBasesInReferenceFull(){
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(hg19MiniReference));
        final CountBasesInReference walker = new CountBasesInReference();
        walker.instanceMain(args.getArgsArray());
        final long[] baseCounts = walker.baseCounts;

        Assert.assertEquals(baseCounts['A'], 4546);
        Assert.assertEquals(baseCounts['C'], 4995);
        Assert.assertEquals(baseCounts['G'], 4559);
        Assert.assertEquals(baseCounts['T'], 3900);
        Assert.assertEquals(baseCounts['N'], 46000);
    }
}