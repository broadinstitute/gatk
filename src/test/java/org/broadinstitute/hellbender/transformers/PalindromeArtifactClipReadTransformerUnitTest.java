package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.stream.IntStream;

import static org.testng.Assert.*;

public class PalindromeArtifactClipReadTransformerUnitTest {

    @Test
    public void test() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final String contig = header.getSequence(0).toString();
        final int refStart = 1000;

        // some random reference sequence -- no ITR yet
        final byte[] refBases = ("GATTCCCCAAGGGGGCTGCTCCCAGAGGGTGTGTTGCTGGGATTGCCCAGGACAGGGATGGCCCTCTCATCAGGTGGGGG" +
                "GCCTCCTCTCGCCGCAGGTCTGGCTGGATGAAGGGCACGGCATAGGTCTGACCTGCCAGGGAGTGCTGCATCCTCACAGG" +
                "GCACAAAAAATGTGCACACACGGGTTCTTCCCACTTTAACCCCTGAGGAATCTGAGGCCTGCTCCTGAAACAGACTGGGC" +
                "GGATGAGTGTGGCATGAAGGGCCTAGGAGATTTCACTTGGGTTTAAAATGCTGTGACCTTGAGTAAGTTGCCGTCTCTGA" +
                "CTGATGCCCAGCCAACTGAGAAACCCAACCCTCTGAGACCAGCACACCCCTTTCAAGCATGTTCCTCCCTCCCCTTCTTT").getBytes();
        final int refEnd = refStart + refBases.length - 1; // inclusive

        final ReferenceDataSource referenceDataSource =
                ReferenceDataSource.of(new ReferenceBases(refBases, new SimpleInterval(contig, refStart, refEnd)), header.getSequenceDictionary());

        final ReadTransformer transformer = new PalindromeArtifactClipReadTransformer(referenceDataSource, Mutect2Engine.MIN_PALINDROME_SIZE);

        // put in an ITR by hand
        final int itrLength = 15;
        // position in reference where first inverted tandem repeat begins and second ends
        final int firstITRStartOffset = 30;
        final int secondITREndOffset = 200;
        for (int n = 0; n < itrLength; n++) {
            refBases[firstITRStartOffset + n] = BaseUtils.simpleComplement(refBases[secondITREndOffset - n]);
        }

        final int numSoftClippedBases = 10;
        final int numUnclippedBases = 30;
        final int readLength = numSoftClippedBases + numUnclippedBases;

        // the fragment spans the ITR region with possibly artifactual soft-clipped bases in the read and the same number of bases in the mate
        final int fragmentLength = (secondITREndOffset - firstITRStartOffset + 1) + 2 * numSoftClippedBases;
        final int readOffset = firstITRStartOffset - numSoftClippedBases;
        final int readStart = readOffset + refStart;

        final byte[] quals = new byte[readLength];
        Arrays.fill(quals, (byte) 30);


        final String cigar = numSoftClippedBases + "S" + numUnclippedBases + "M";

        final byte[] readBasesMatchingRef = Arrays.copyOfRange(refBases, readOffset, readOffset + readLength);
        final byte[] readBasesWithArtifact = Arrays.copyOf(readBasesMatchingRef, readLength);
        for (int n = 0; n < numSoftClippedBases; n++) {
            readBasesWithArtifact[n] = BaseUtils.simpleComplement(refBases[readOffset + fragmentLength - 1 - n]);
        }


        // read is upstream of mate and has the artifact
        final GATKRead artifactRead = makeRead(header, contig, readStart, fragmentLength, readBasesWithArtifact, quals, cigar);
        Assert.assertTrue(transformer.apply(artifactRead).getCigar().getFirstCigarElement().getOperator() == CigarOperator.HARD_CLIP);

        // shift the mate by one base and no hard-clipping should occur
        final GATKRead notQuiteArtifactRead = makeRead(header, contig, readStart, fragmentLength + 1, readBasesWithArtifact, quals, cigar);
        Assert.assertTrue(transformer.apply(notQuiteArtifactRead).getCigar().getFirstCigarElement().getOperator() == CigarOperator.SOFT_CLIP);

        // definitely no artifact
        final GATKRead definitelyNotArtifactRead = makeRead(header, contig, readStart, fragmentLength, readBasesMatchingRef, quals, cigar);
        Assert.assertTrue(transformer.apply(definitelyNotArtifactRead).getCigar().getFirstCigarElement().getOperator() == CigarOperator.SOFT_CLIP);

    }

    private static GATKRead makeRead(final SAMFileHeader header, final String contig, final int readStart, final int fragmentLength, final byte[] bases, final byte[] qual, final String cigar) {
        final SAMRecord read = new SAMRecord(header);
        read.setReferenceName(contig);
        read.setAlignmentStart(readStart);
        read.setReadPairedFlag(true);
        read.setReadUnmappedFlag(false);
        read.setMateUnmappedFlag(false);
        read.setMateReferenceName("Mate");
        read.setReadNegativeStrandFlag(false);
        read.setMateNegativeStrandFlag(true);
        read.setReadBases(bases);
        read.setBaseQualities(qual);


        final int mateStart = readStart + fragmentLength - bases.length;
        read.setMateAlignmentStart(mateStart);
        read.setInferredInsertSize(fragmentLength);
        read.setProperPairFlag(true);
        read.setCigarString(cigar);


        return new SAMRecordToGATKReadAdapter(read);
    }

}