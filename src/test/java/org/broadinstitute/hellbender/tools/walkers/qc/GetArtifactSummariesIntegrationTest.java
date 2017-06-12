package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.junit.Test;
import org.testng.Assert;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

/**
 * Created by David Benjamin on 6/14/17.
 */
public class GetArtifactSummariesIntegrationTest extends CommandLineProgramTest  {

    @Test
    public void test() {
        final File NA12878 = new File(largeFileTestDir, "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam");

        final File output = createTempFile("output", ".table");
        final String[] args = {
                "-I", NA12878.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "-L", "20:9999950-10000125",
                "-O", output.getAbsolutePath(),
                "-" + GetArtifactSummaries.BASES_BEFORE_AND_AFTER_SHORT_NAME, "3",
                "-" + GetArtifactSummaries.MIN_BASE_QUALITY_SHORT_NAME, "20"
        };
        runCommandLine(args);

        final List<ArtifactSummary> result = ArtifactSummary.readArtifactSummaries(output);

        /*
        true counts determined manually in IGV
        9,999,976:  a low-quality T->C and a soft clip.  Shouldn't call anything.
        9,999,990: two low-quality substitutions.  Shouldn't call anything.
        9,999,996: everything is insertion start
        10,000,117: het C/T
         */
        final int[] positions = new int[] { 9999976, 9999990, 9999996, 10000117};
        final String[] expectedRefBases = new String[] {"CTCTATC", "CAGTACC", "CACAGTT", "GGACTTT" };
        final int[][] expectedACGTCounts = new int[][] {{0,0,0,68}, {0,0,0,81}, {0,0,0,0}, {0,35,0,28} };
        final int[] expectedInsertionStartCounts = new int[] {0, 0, 83, 0 };
        final int[] expectedDeletionStartCounts = new int[] {0, 0, 0, 0 };
        final int[] indicesToVerify = Arrays.stream(positions).map(n -> n - 9999950).toArray();

        for (int n : new int[] {0,1,2,3}) {
            final ArtifactSummary artifactSummary = result.get(indicesToVerify[n]);
            Assert.assertEquals(artifactSummary.getContig(), "20");
            Assert.assertEquals(artifactSummary.getRefBases(), expectedRefBases[n].getBytes());
            Assert.assertEquals(artifactSummary.getStart(), positions[n]);
            Assert.assertEquals(artifactSummary.getBaseCounts(), expectedACGTCounts[n]);
            Assert.assertEquals(artifactSummary.getDeletionStartCount(), expectedDeletionStartCounts[n]);
            Assert.assertEquals(artifactSummary.getInsertionStartCount(), expectedInsertionStartCounts[n]);
        }
    }
}