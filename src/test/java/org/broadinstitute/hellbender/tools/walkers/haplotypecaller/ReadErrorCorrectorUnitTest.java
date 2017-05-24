package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class ReadErrorCorrectorUnitTest {
    private static final boolean debug = true;
    final String refChunk = "GCATAAACATGGCTCACTGC";
    final String refChunkHard = "AGCCTTGAACTCCTGGGCTCAAGTGATCCTCCTGCCTCAGTTTCCCATGTAGCTGGGACCACAGGTGGGGGCTCCACCCCTGGCTGATTTTTTTTTTTTTTTTTTTTTGAGATAGGGT";

    @Test
    public void TestBasicCorrectionSet() {

        final byte[] trueBases = refChunk.getBytes();
        final int numCorrections = 50;
        final ReadErrorCorrector.CorrectionSet correctionSet = new ReadErrorCorrector.CorrectionSet(trueBases.length);

        int offset = 2;
        for (int k=0; k < numCorrections; k++) {
            // introduce one correction at a random offset in array. To make testing easier, we will replicate corrrection
            final byte base = trueBases[offset];
            correctionSet.add(offset, base);
            // skip to some other offset
            offset += 7;
            if (offset >= trueBases.length)
                offset -= trueBases.length;
        }

        for (int k=0; k < trueBases.length; k++) {
            final byte corr =  correctionSet.getConsensusCorrection(k);
            Assert.assertEquals(corr, trueBases[k]);
        }
    }

    @Test
    public void TestExtendedCorrectionSet() {

        final byte[] trueBases = refChunk.getBytes();
        final int numCorrections = 50;
        final ReadErrorCorrector.CorrectionSet correctionSet = new ReadErrorCorrector.CorrectionSet(trueBases.length);

        for (int offset=0; offset < trueBases.length; offset++) {
            // insert k corrections at offset k and make sure we get exactly k bases back
            for (int k=0; k < offset; k++)
                correctionSet.add(offset,trueBases[offset]);

        }

        for (int offset=0; offset < trueBases.length; offset++) {
            Assert.assertEquals(correctionSet.get(offset).size(),offset);
        }
    }

    @Test
    public void TestAddReadsToKmers() {
        final int NUM_GOOD_READS = 500;

        final String bases = "AAAAAAAAAAAAAAA";
        final int READ_LENGTH = bases.length();
        final int kmerLengthForReadErrorCorrection = READ_LENGTH;
        final List<GATKRead> finalizedReadList = new ArrayList<GATKRead>(NUM_GOOD_READS);
        int offset = 0;
        final byte[] quals = new byte[READ_LENGTH];

        Arrays.fill(quals,(byte)30);

        for (int k=0; k < NUM_GOOD_READS; k++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(bases.getBytes(), quals,READ_LENGTH+"M");
            finalizedReadList.add(read);
        }

        ReadErrorCorrector readErrorCorrector = new ReadErrorCorrector(kmerLengthForReadErrorCorrection,(byte)6,10, debug,refChunkHard.getBytes());
        readErrorCorrector.addReadsToKmers(finalizedReadList);

        // special trivial case: kmer length is equal to read length.
        // K-mer counter should hold then exactly one kmer
        Assert.assertEquals(readErrorCorrector.countsByKMer.getCountedKmers().size(), 1);
        for (final KMerCounter.CountedKmer kmer : readErrorCorrector.countsByKMer.getCountedKmers()) {
            Assert.assertTrue(Arrays.equals( kmer.getKmer().bases(),bases.getBytes()));
            Assert.assertEquals(kmer.getCount(),NUM_GOOD_READS);
        }

        // special case 2: kmers are all the same but length < read length.
        // Each kmer is added then readLength-kmerLength+1 times
        final int KMER_LENGTH = 10;
        readErrorCorrector = new ReadErrorCorrector(KMER_LENGTH,(byte)6,10, debug,refChunkHard.getBytes());
        readErrorCorrector.addReadsToKmers(finalizedReadList);
        Assert.assertEquals(readErrorCorrector.countsByKMer.getCountedKmers().size(), 1);
        for (final KMerCounter.CountedKmer kmer : readErrorCorrector.countsByKMer.getCountedKmers()) {
            Assert.assertEquals(kmer.getCount(),NUM_GOOD_READS*(READ_LENGTH-KMER_LENGTH+1));
        }

    }
    @Test
    public void TestBasicErrorCorrection() {
        final int NUM_GOOD_READS = 500;
        final int NUM_BAD_READS = 10;
        final int READ_LENGTH = 15;
        final int kmerLengthForReadErrorCorrection = 10;
        final List<GATKRead> finalizedReadList = new ArrayList<GATKRead>(NUM_GOOD_READS);
        int offset = 0;
        final byte[] quals = new byte[READ_LENGTH];

        Arrays.fill(quals,(byte)30);

        for (int k=0; k < NUM_GOOD_READS; k++) {
            final byte[] bases = Arrays.copyOfRange(refChunk.getBytes(),offset,offset+READ_LENGTH);
            final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals,READ_LENGTH+"M");
            finalizedReadList.add(read);
            offset++;
            if (offset >= refChunk.length()-READ_LENGTH)
                offset = 0;
        }
        offset = 2;
        // coverage profile is now perfectly triangular with "good" bases. Inject now bad bases with errors in them.
        for (int k=0; k < NUM_BAD_READS; k++) {
            final byte[] bases = finalizedReadList.get(k).getBases().clone();
            bases[offset] = 'N';
            final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, READ_LENGTH + "M");
            finalizedReadList.add(read);
            offset += 7;
            if (offset >= READ_LENGTH)
                offset = 4; // just some randomly circulating offset for error position
        }

        // now correct all reads
        final ReadErrorCorrector readErrorCorrector = new ReadErrorCorrector(kmerLengthForReadErrorCorrection,(byte)6,10, debug,refChunkHard.getBytes());
        readErrorCorrector.addReadsToKmers(finalizedReadList);
        readErrorCorrector.correctReads(finalizedReadList);

        // check that corrected reads have exactly same content as original reads
        for (int k=0; k < NUM_BAD_READS; k++) {
            final byte[] badBases = finalizedReadList.get(k).getBases();
            final byte[] originalBases = finalizedReadList.get(k).getBases();
            Assert.assertTrue(Arrays.equals(badBases,originalBases));
        }
    }
}
