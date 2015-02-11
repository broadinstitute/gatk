package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.sam.ArtificialSAMUtils;

import java.util.List;

public class LIBSTest {
    public static final int locus = 44367788;
    final String cigarString;
    final int readLength;
    final private List<CigarElement> elements;

    public LIBSTest(final String cigarString) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        this.cigarString = cigarString;
        this.elements = cigar.getCigarElements();
        this.readLength = cigar.getReadLength();
    }

    @Override
    public String toString() {
        return "LIBSTest{" +
                "cigar='" + cigarString + '\'' +
                ", readLength=" + readLength +
                '}';
    }

    public List<CigarElement> getElements() {
        return elements;
    }

    public SAMRecord makeRead() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, locus, readLength);
        read.setReadBases(Utils.dupBytes((byte) 'A', readLength));
        final byte[] quals = new byte[readLength];
        for ( int i = 0; i < readLength; i++ )
            quals[i] = (byte)(i % QualityUtils.MAX_SAM_QUAL_SCORE);
        read.setBaseQualities(quals);
        read.setCigarString(cigarString);
        return read;
    }
}
