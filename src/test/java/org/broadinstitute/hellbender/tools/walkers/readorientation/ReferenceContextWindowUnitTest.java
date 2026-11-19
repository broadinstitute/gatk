package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.Random;

public class ReferenceContextWindowUnitTest {
    private static final String CONTIG = "1";

    @DataProvider
    public Object[][] windowSizes() {
        return new Object[][] {{3}, {4}, {7}, {100}, {ReferenceContextWindow.DEFAULT_WINDOW_SIZE}};
    }

    @Test(dataProvider = "windowSizes")
    public void testContextsAgreeWithGetKmerAround(final int windowSize) {
        // a contig with Ns at both ends, in runs and alone, so that every position's context, including those running
        // off the ends of the contig, can be compared with what F1R2CountsCollector used to take from getKmerAround
        final Random random = new Random(windowSize);
        final byte[] bases = new byte[2_500];
        for (int n = 0; n < bases.length; n++) {
            bases[n] = (byte) "ACGTACGTACGTACGTN".charAt(random.nextInt(17));
        }
        bases[0] = 'N';
        bases[bases.length - 1] = 'N';
        for (int n = 1_000; n < 1_020; n++) {
            bases[n] = 'N';
        }
        final ReferenceMemorySource source = new ReferenceMemorySource(new ReferenceBases(bases, new SimpleInterval(CONTIG, 1, bases.length)),
                new SAMSequenceDictionary(Collections.singletonList(new SAMSequenceRecord(CONTIG, bases.length))));

        final ReferenceContextWindow window = new ReferenceContextWindow(windowSize);
        int numContexts = 0;
        for (int position = 1; position <= bases.length; position++) {
            final ReferenceContext referenceContext = new ReferenceContext(source, new SimpleInterval(CONTIG, position, position));
            final String kmer = referenceContext.getKmerAround(position, F1R2FilterConstants.REF_CONTEXT_PADDING);
            final String expected = kmer == null || kmer.contains("N") ? null : kmer;
            Assert.assertEquals(window.contextAt(referenceContext, CONTIG, position), expected, "position " + position);
            numContexts += expected == null ? 0 : 1;
        }
        Assert.assertTrue(numContexts > 1_000);
    }
}
