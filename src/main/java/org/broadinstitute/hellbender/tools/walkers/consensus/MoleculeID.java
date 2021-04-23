package org.broadinstitute.hellbender.tools.walkers.consensus;

import htsjdk.samtools.SAMTag;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * A container class for the molecule ID, which consists of an integer ID and a binary strand.
 * For example, Reads with the tags 12/A and 12/B originated from the same DNA fragment before PCR,
 * (i.e. from the same library) but they originated from different strands in that library.
 * In other words, one read is F1R2 and the other F2R1.
 *
 * The word "molecule" here refers to the original DNA fragment with barcode before undergoing
 * PCR and sequencing. We amplify this molecule through PCR and end up with many duplicate _fragments_.
 */
public class MoleculeID {
    private int moleculeNumber;
    private String strand;

    public MoleculeID(final GATKRead read){
        this.moleculeNumber = getMoleculeNumberOfRead(read);
        this.strand = getStrandOfRead(read);
    }

    public MoleculeID(final int moleculeNumber, final String strand){
        this.moleculeNumber = moleculeNumber;
        this.strand = strand;
    }

    public int getMoleculeNumber() {
        return moleculeNumber;
    }

    public String getStrand() {
        return strand;
    }

    /** Format the molecule ID as stored in the sam/bam/cram file under the {@link SAMTag.MI} tag **/
    public String getSAMField(){
        return moleculeNumber + ReadsWithSameUMI.FGBIO_MI_TAG_DELIMITER + strand;
    }

    /** Extracts the molecule number portion of the {@link SAMTag.MI} field of the read **/
    public static int getMoleculeNumberOfRead(final GATKRead read){
        final String MITag = read.getAttributeAsString(SAMTag.MI.name());
        return Integer.parseInt(MITag.split(ReadsWithSameUMI.FGBIO_MI_TAG_DELIMITER)[0]);
    }

    /** Extracts the strand portion of the {@link SAMTag.MI} field of the read **/
    public static String getStrandOfRead(final GATKRead read){
        final String MITag = read.getAttributeAsString(SAMTag.MI.name());
        return MITag.split(ReadsWithSameUMI.FGBIO_MI_TAG_DELIMITER)[1];
    }

    /**
     * Assumes that the input reads have the same molecule number in the {@link SAMTag.MI} tag
     * @returns Counts of reads from each strand, the first element is always larger than the second
     **/
    public static Pair<Integer, Integer> countStrands(final List<GATKRead> reads){
        final int strandACount = (int) reads.stream().filter(r -> getStrandOfRead(r).equals("A")).count();
        final int strandBCount = (int) reads.stream().filter(r -> getStrandOfRead(r).equals("B")).count();
        return new ImmutablePair<>(strandACount, strandBCount);
    }
}
