package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Map;

public class ExampleReferenceWalkerIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testExampleReferenceWalker(){
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference( getTestFile(("example.fasta")));
        args.addInput(getTestFile("veryspecific.bam"));
        args.addVCF(getTestFile("example.vcf"));

        final ExampleReferenceWalker exampleReferenceWalker = new ExampleReferenceWalker();
        exampleReferenceWalker.instanceMain(args.getArgsArray());

        final Map<String, ExampleReferenceWalker.OverlapCounts> contextCounts = exampleReferenceWalker.contextCounts;
        final ExampleReferenceWalker.OverlapCounts aaaCounts = contextCounts.get("AAA");
        final ExampleReferenceWalker.OverlapCounts nnnCounts = contextCounts.get("NNN");
        final ExampleReferenceWalker.OverlapCounts aacCounts = contextCounts.get("AAC");
        assertCountsEqual(aaaCounts, 78, 9, 3 );
        assertCountsEqual(nnnCounts, 78, 1, 0);
        assertCountsEqual(aacCounts, 1, 0, 1);
    }

    private static void assertCountsEqual(ExampleReferenceWalker.OverlapCounts counts, long seen, long reads, long variants){
        Assert.assertEquals(counts.timesSeen, seen);
        Assert.assertEquals(counts.overlappedByReads, reads);
        Assert.assertEquals(counts.overlappedByVariants, variants);
    }
}
