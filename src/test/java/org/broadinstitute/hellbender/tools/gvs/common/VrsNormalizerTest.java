package org.broadinstitute.hellbender.tools.gvs.common;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.Iterator;

public class VrsNormalizerTest extends GATKBaseTest {

    private static VrsNormalizer normalizer;
    private static ReferenceContext testRefContext;

    private static class TestReferenceDataSource implements ReferenceDataSource {
        private final String chr20Sequence;

        public TestReferenceDataSource() {
            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < 1000; i++) {
                builder.append('C');
            }

            // For substitution trimming test: [200,204) = ACGT
            builder.setCharAt(200, 'A');
            builder.setCharAt(201, 'C');
            builder.setCharAt(202, 'G');
            builder.setCharAt(203, 'T');

            // For no-variation test: [300,302) = AC
            builder.setCharAt(300, 'A');
            builder.setCharAt(301, 'C');

            // For deletion test: [400,402) = CA
            builder.setCharAt(400, 'C');
            builder.setCharAt(401, 'A');

            // For unambiguous insertion test around index 500 (prev != T, current != A)
            builder.setCharAt(499, 'C');
            builder.setCharAt(500, 'G');

            // For reference-derived ambiguous insertion test around boundary index 601:
            // ... T A T A T A ...
            builder.setCharAt(600, 'T');
            builder.setCharAt(601, 'A');
            builder.setCharAt(602, 'T');
            builder.setCharAt(603, 'A');
            builder.setCharAt(604, 'T');
            builder.setCharAt(605, 'A');

            chr20Sequence = builder.toString();
        }

        @Override
        public ReferenceSequence queryAndPrefetch(String contig, long start, long stop) {
            if (!"chr20".equals(contig)) {
                throw new IllegalArgumentException("Unsupported contig in test data source: " + contig);
            }
            int zeroBasedStart = (int) start - 1;
            int zeroBasedStopExclusive = (int) stop;
            if (zeroBasedStart < 0 || zeroBasedStopExclusive > chr20Sequence.length() || zeroBasedStart > zeroBasedStopExclusive) {
                throw new IllegalArgumentException("Invalid query range for test sequence: " + start + "-" + stop);
            }
            String seq = chr20Sequence.substring(zeroBasedStart, zeroBasedStopExclusive);
            return new ReferenceSequence(contig, (int) start, seq.getBytes());
        }

        @Override
        public ReferenceSequence queryAndPrefetch(SimpleInterval interval) {
            return queryAndPrefetch(interval.getContig(), interval.getStart(), interval.getEnd());
        }

        @Override
        public SAMSequenceDictionary getSequenceDictionary() {
            SAMSequenceDictionary dict = new SAMSequenceDictionary();
            dict.addSequence(new SAMSequenceRecord("chr20", 63025520));
            return dict;
        }

        @Override
        public void close() {
            // no-op
        }

        @Override
        public Iterator<Byte> iterator() {
            throw new UnsupportedOperationException("iterator() not implemented for test");
        }
    }

    @BeforeClass
    public void setUp() {
        TestReferenceDataSource refSource = new TestReferenceDataSource();
        normalizer = new VrsNormalizer();
        // Wrap the test ReferenceDataSource in a ReferenceContext centred on chr20
        // with a large window so all test positions are covered.
        testRefContext = new ReferenceContext(refSource, new SimpleInterval("chr20", 1, 1000), 0, 0);
    }

    /** Helper so tests don't need to repeat the refContext argument. */
    private NormalizedAllele normalize(String contig, long start, long end, String alt) {
        return normalizer.normalize(contig, start, end, alt, testRefContext);
    }

    @Test
    public void testSnpPassThroughAsLiteralSequence() {
        NormalizedAllele allele = normalize("chr20", 99L, 100L, "T");
        Assert.assertNotNull(allele);
        Assert.assertEquals(allele.start, 99L);
        Assert.assertEquals(allele.end, 100L);
        Assert.assertEquals(allele.state.type, NormalizedAllele.SequenceState.StateType.LITERAL_SEQUENCE_EXPRESSION);
        Assert.assertEquals(allele.state.sequence, "T");
        Assert.assertNull(allele.state.length);
        Assert.assertNull(allele.state.repeatSubunitLength);
    }

    @Test
    public void testSubstitutionWithTrimmingReturnsTrimmedCoordinatesAndState() {
        NormalizedAllele allele = normalize("chr20", 200L, 204L, "ATGT");

        Assert.assertNotNull(allele);
        Assert.assertEquals(allele.start, 201L);
        Assert.assertEquals(allele.end, 202L);
        Assert.assertEquals(allele.state.type, NormalizedAllele.SequenceState.StateType.LITERAL_SEQUENCE_EXPRESSION);
        Assert.assertEquals(allele.state.sequence, "T");
    }

    @Test
    public void testNoVariationReturnsReferenceLengthState() {
        NormalizedAllele allele = normalize("chr20", 300L, 302L, "AC");

        Assert.assertNotNull(allele);
        Assert.assertEquals(allele.start, 300L);
        Assert.assertEquals(allele.end, 302L);
        Assert.assertEquals(allele.state.type, NormalizedAllele.SequenceState.StateType.REFERENCE_LENGTH_EXPRESSION);
        Assert.assertEquals(allele.state.length, Integer.valueOf(2));
        Assert.assertEquals(allele.state.repeatSubunitLength, Integer.valueOf(2));
    }

    @Test
    public void testDeletionReturnsReferenceLengthState() {
        NormalizedAllele allele = normalize("chr20", 400L, 402L, "C");

        Assert.assertNotNull(allele);
        Assert.assertEquals(allele.state.type, NormalizedAllele.SequenceState.StateType.REFERENCE_LENGTH_EXPRESSION);
        Assert.assertEquals(allele.state.length, Integer.valueOf(1));
        Assert.assertEquals(allele.state.repeatSubunitLength, Integer.valueOf(1));
    }

    @Test
    public void testUnambiguousInsertionReturnsLiteralState() {
        NormalizedAllele allele = normalize("chr20", 500L, 500L, "AT");

        Assert.assertNotNull(allele);
        Assert.assertEquals(allele.start, 500L);
        Assert.assertEquals(allele.end, 500L);
        Assert.assertEquals(allele.state.type, NormalizedAllele.SequenceState.StateType.LITERAL_SEQUENCE_EXPRESSION);
        Assert.assertEquals(allele.state.sequence, "AT");
    }

    @Test
    public void testAmbiguousInsertionReferenceDerivedReturnsReferenceLengthState() {
        NormalizedAllele allele = normalize("chr20", 601L, 601L, "AT");

        Assert.assertNotNull(allele);
        Assert.assertEquals(allele.start, 600L);
        Assert.assertEquals(allele.end, 606L);
        Assert.assertEquals(allele.state.type, NormalizedAllele.SequenceState.StateType.REFERENCE_LENGTH_EXPRESSION);
        Assert.assertEquals(allele.state.length, Integer.valueOf(8));
        Assert.assertEquals(allele.state.repeatSubunitLength, Integer.valueOf(2));
    }

    @Test
    public void testAmbiguousInsertionNovelReturnsLiteralState() {
        NormalizedAllele allele = normalize("chr20", 601L, 601L, "AG");

        Assert.assertNotNull(allele);
        Assert.assertEquals(allele.start, 601L);
        Assert.assertEquals(allele.end, 602L);
        Assert.assertEquals(allele.state.type, NormalizedAllele.SequenceState.StateType.LITERAL_SEQUENCE_EXPRESSION);
        Assert.assertEquals(allele.state.sequence, "AGA");
    }
}