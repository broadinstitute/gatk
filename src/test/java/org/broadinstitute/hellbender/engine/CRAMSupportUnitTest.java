package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.cram.build.CRAMReferenceRegion;
import htsjdk.samtools.cram.ref.ReferenceContext;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.cram.structure.AlignmentContext;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CRAMSupportUnitTest extends GATKBaseTest {

    private static final int REFERENCE_SEQUENCE_ZERO = 0;
    private static final int REFERENCE_CONTIG_LENGTH = 10000;

    private static final SAMFileHeader SAM_FILE_HEADER = createSAMFileHeader();

    private static SAMFileHeader createSAMFileHeader() {
        final List<SAMSequenceRecord> sequenceRecords = new ArrayList<>();
        sequenceRecords.add(new SAMSequenceRecord("0", REFERENCE_CONTIG_LENGTH));
        sequenceRecords.add(new SAMSequenceRecord("1", REFERENCE_CONTIG_LENGTH));
        final SAMFileHeader header = new SAMFileHeader(new SAMSequenceDictionary(sequenceRecords));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    /*
     * This is a unit test ported from CRAMReferenceRegionTest in HTSJDK that tests for the conditions that
     * trigger the CRAM base corruption bug reported in https://github.com/broadinstitute/gatk/issues/8768 and
     * fixed in HTSJDK 4.1.1 / GATK 4.6.0.0.
     *
     * This test fails using HTSJDK version 4.1.0, which predates the fix in https://github.com/samtools/htsjdk/pull/1708
     * for the bug reported in https://github.com/broadinstitute/gatk/issues/8768
     *
     * Simulates the state transitions that occur when writing a CRAM file; specifically, use transitions that mirror
     * the ones that occur when writing a CRAM with the conditions that affect (and that fail without the fix to)
     * https://github.com/broadinstitute/gatk/issues/8768, i.e., a container with one or more containers with reads
     * aligned to position 1 of a given contig, followed by one or more containers with reads aligned past position 1
     * of the same contig.
     */
    @Test
    public void testCRAMBaseCorruptionBugIssue8768() {
        // start with an entire reference sequence
        final CRAMReferenceRegion cramReferenceRegion = getAlternatingReferenceRegion();
        cramReferenceRegion.fetchReferenceBases(REFERENCE_SEQUENCE_ZERO);
        final long fullRegionFragmentLength = cramReferenceRegion.getRegionLength();
        Assert.assertEquals(fullRegionFragmentLength, REFERENCE_CONTIG_LENGTH);

        // transition to a shorter reference fragment using fetchReferenceBasesByRegion, then back to the full region
        final int SHORT_FRAGMENT_LENGTH = 5;
        Assert.assertTrue(SHORT_FRAGMENT_LENGTH < fullRegionFragmentLength);
        cramReferenceRegion.fetchReferenceBasesByRegion(REFERENCE_SEQUENCE_ZERO, 0, SHORT_FRAGMENT_LENGTH);
        Assert.assertEquals(cramReferenceRegion.getRegionLength(), SHORT_FRAGMENT_LENGTH);

        // now transition back to the full sequence; this is where the bug previously would have occurred
        cramReferenceRegion.fetchReferenceBases(REFERENCE_SEQUENCE_ZERO);
        // this assert would fail without the fix
        Assert.assertEquals(cramReferenceRegion.getRegionLength(), fullRegionFragmentLength);

        // transition to a shorter region fragment length using fetchReferenceBasesByRegion(AlignmentContext), then back to the full region
        Assert.assertTrue(SHORT_FRAGMENT_LENGTH < fullRegionFragmentLength);
        cramReferenceRegion.fetchReferenceBasesByRegion(
                new AlignmentContext(
                        new ReferenceContext(REFERENCE_SEQUENCE_ZERO),
                        1,
                        SHORT_FRAGMENT_LENGTH));
        Assert.assertEquals(cramReferenceRegion.getRegionLength(), SHORT_FRAGMENT_LENGTH);

        // now transition back to the full sequence
        cramReferenceRegion.fetchReferenceBases(REFERENCE_SEQUENCE_ZERO);
        Assert.assertEquals(cramReferenceRegion.getRegionLength(), fullRegionFragmentLength);
    }

    private static CRAMReferenceRegion getAlternatingReferenceRegion() {
        return new CRAMReferenceRegion(
                new ReferenceSource(getReferenceFileWithAlternatingBases(REFERENCE_CONTIG_LENGTH)),
                SAM_FILE_HEADER.getSequenceDictionary());
    }

    private static InMemoryReferenceSequenceFile getReferenceFileWithAlternatingBases(final int length) {
        final InMemoryReferenceSequenceFile referenceFile = new InMemoryReferenceSequenceFile();

        // one contig with repeated ACGT...
        final byte[] seq0Bases = getRepeatingBaseSequence(REFERENCE_CONTIG_LENGTH, false);
        referenceFile.add("0", seq0Bases);

        // one contig with repeated TGCA...
        final byte[] seq1Bases = getRepeatingBaseSequence(REFERENCE_CONTIG_LENGTH, true);
        referenceFile.add("1", seq1Bases);

        return referenceFile;
    }

    // fill an array with the repeated base sequence "ACGTACGTACGT...", or reversed
    private static byte[] getRepeatingBaseSequence(final int length, final boolean reversed) {
        byte[] bases = new byte[length];
        for (int i = 0; (i + 4) < bases.length; i += 4) {
            bases[i] = (byte) (reversed ? 'T' : 'A');
            bases[i+1] = (byte) (reversed ? 'G' : 'C');
            bases[i+2] = (byte) (reversed ? 'C' : 'G');
            bases[i+3] = (byte) (reversed ? 'A' : 'T');
        }
        return bases;
    }
    private static class InMemoryReferenceSequenceFile implements
            ReferenceSequenceFile {
        Map<String, ReferenceSequence> map = new HashMap<String, ReferenceSequence>();
        List<String> index;
        int current = 0;

        public void add(final String name, final byte[] bases) {
            final ReferenceSequence sequence = new ReferenceSequence(name,
                    map.size(), bases);
            map.put(sequence.getName(), sequence);
        }

        @Override
        public void reset() {
            current = 0;
        }

        @Override
        public ReferenceSequence nextSequence() {
            if (current >= index.size()) return null;
            return map.get(index.get(current++));
        }

        @Override
        public boolean isIndexed() {
            return true;
        }

        @Override
        public ReferenceSequence getSubsequenceAt(final String contig, final long start,
                                                  final long stop) {
            final ReferenceSequence sequence = getSequence(contig);
            if (sequence == null) return null;
            final byte[] bases = new byte[(int) (stop - start) + 1];
            System.arraycopy(sequence.getBases(), (int) start - 1, bases, 0,
                    bases.length);
            return new ReferenceSequence(contig, sequence.getContigIndex(),
                    bases);
        }

        @Override
        public SAMSequenceDictionary getSequenceDictionary() {
            return null;
        }

        @Override
        public ReferenceSequence getSequence(final String contig) {
            return map.get(contig);
        }

        @Override
        public void close() throws IOException {
            map.clear();
            index.clear();
            current = 0;
        }
    }
}
