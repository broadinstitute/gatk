package org.broadinstitute.hellbender.tools.walkers.consensus;

import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A container class for a set of reads that share the same unique molecular identifier (UMI) as judged by
 * FGBio GroupReadsByUmi (http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html)
 *
 * Examples of molecule IDs (MI tag):
 *
 * "0/A" (The first molecule in the bam, A strand)
 * "0/B" (The first molecule in the bam, B strand)
 * "99/A" (100th molecule in the bam, A strand)
 *
 * For a given set of reads with the same molecule number, the strand with a larger number of reads is defined as the A strand.
 * i.e. A and B strand doesn't map to top or bottom strand.
 *
 * I use "top strand" as a synonym for "F1R2". "Bottom strand" is synonymous to "F2R1."
 *
 * Thus only the integer component is relevant for identifying reads that originated from the same molecule.
 * Should the need arise, we could extend this to distinguish between different strands of the same molecule.
 *
 * All reads in the set must share the same molecular number; this allows for reads that originated from the same fragment before PCR
 * but from the different strands to be grouped in the same set.
 * For instance, 1/A and 1/B may be in the same set, as they share the same UMI.
 * But 2/A and 3/A may not be in the same set.
 */
public class ReadsWithSameUMI implements Locatable {
    public static final String FGBIO_MI_TAG_DELIMITER = "/";
    public MoleculeID moleculeID;

    private SimpleInterval interval;

    private List<GATKRead> reads;

    public ReadsWithSameUMI(final GATKRead read){
        reads = new ArrayList<>();
        init(read);
    }

    private void init(GATKRead read){
        Utils.validate(reads.isEmpty(), String.format("Initializing a non-empty set"));
        moleculeID = new MoleculeID(read);
        interval = new SimpleInterval(read);
        reads.add(read);
    }

    public List<GATKRead> getReads(){
        return Collections.unmodifiableList(reads);
    }

    /**
     * Add a read to the set. Throws an error if the molecule ID doens't match.
     * **/
    public void addRead(final GATKRead read){
        Utils.validate(reads.isEmpty() || moleculeID.getMoleculeNumber() == MoleculeID.getMoleculeNumberOfRead(read),
                String.format("Molecule number of the set and that of the new read don't match: set number = %d, read number = %d", moleculeID.getMoleculeNumber(), MoleculeID.getMoleculeNumberOfRead(read)));
        Utils.validate(interval.contigsMatch(read),
                String.format("Adding a read from another contig: set contig = %s, read contig = %s", interval.getContig(), read.getContig()));

        if (reads.isEmpty()){
            init(read);
            return;
        }

        reads.add(read);
        interval = interval.spanWith(read);
    }

    public boolean isEmpty(){
        return reads.isEmpty();
    }

    public int getMoleculeNumber() {
        Utils.nonNull(moleculeID, "Querying a non-initialized set for a molecule number");
        return moleculeID.getMoleculeNumber();
    }

    public SimpleInterval getInterval(){
        return interval;
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }
}
