package org.broadinstitute.hellbender.tools.walkers.mutect.consensus;

import org.broadinstitute.hellbender.tools.walkers.mutect.UMI;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class DuplicateSet {
    public static final String FGBIO_MOLECULAR_IDENTIFIER_TAG = "MI";
    public static final String FGBIO_MI_TAG_DELIMITER = "/";
    private int moleculeId = -1; // TODO: extract an ID class.
    private UMI umi;
    private String contig;
    private int fragmentStart = -1;
    private int fragmentEnd = -1;
    private List<GATKRead> reads;
    boolean smallInsert; // if true, the reads read into adaptors
    private boolean paired;  // TODO: is this useful/how can I detect this?
    // IN FACT, start here, also write tests --- count the number of MI=1's, MI=2's, etc...
    public DuplicateSet(){
        reads = new ArrayList<>();
    }

    public DuplicateSet(final GATKRead read){
        reads = new ArrayList<>();
        init(read);
        reads.add(read);
    }

    public void init(GATKRead read){
        Utils.validate(moleculeId == -1 || moleculeId == getMoleculeID(read),
                String.format("Inconsisntent molecule IDs: Duplicate set id = %s, read molecule id = %s", moleculeId, getMoleculeID(read)));
        setMoleduleId(read);umi = new UMI(read);
        contig = read.getContig();
        fragmentStart = read.getStart();
        fragmentEnd = read.getEnd(); // TODO: does this include softclips?
        paired = false;
    }

    public List<GATKRead> getReads(){
        return reads;
    }

    public boolean sameMolecule(final GATKRead read){
        return getMoleculeID(read) == moleculeId;
    }

    public void setMoleduleId(GATKRead read){
        moleculeId = getMoleculeID(read);
    }

    public void setMoleduleId(int id){
        moleculeId = id;
    }

    /**
     * Some examples of molecule IDs (MI tag):
     *
     * "0/A" (The first molecule in the bam, top (A) strand)
     * "0/B" (The first molecule in the bam, bottom (B) strand)
     * "99/A" (100th molecule in the bam, top top (A) strand)
     *
     * Top strand is synonymous to "F1R2"
     * Bottom strand is synonymous to "F2R1"
     *
     * Thus only the integer component is relevant for identifying reads that originated from the same molecule.
     * Should the need arise, we could extend this to distinguish between different strands of the same molecule.
     */
    public static int getMoleculeID(final GATKRead read) {
        final String MITag = read.getAttributeAsString(FGBIO_MOLECULAR_IDENTIFIER_TAG);
        return Integer.parseInt(MITag.split(FGBIO_MI_TAG_DELIMITER)[0]);
    }

    /** Returns true if the read was properly added to the duplicate set **/
    public boolean addRead(final GATKRead read){
        if (reads.isEmpty()){
            init(read);
            reads.add(read);
            return true;
        }

        if (sameMolecule(read)){
            reads.add(read);
            if (read.getStart() < fragmentStart){
                fragmentStart = read.getStart();
            }
            if (read.getEnd() > fragmentEnd){
                fragmentEnd = read.getEnd();
            }
            return true;
        } else {
            return false;
        }
    }

    public int getFragmentStart(){
        return fragmentStart;
    }

    public int getFragmentEnd(){
        return fragmentEnd;
    }

    public int getMoleculeId() { return moleculeId; }

    public SimpleInterval getDuplicateSetInterval(){
        return new SimpleInterval(contig, fragmentStart, fragmentEnd);
    }

    public boolean hasValidInterval(){
        return SimpleInterval.isValid(contig, fragmentStart, fragmentEnd);
    }

    public static List<String> getMolecularIDs(final List<GATKRead> reads) {
        return reads.stream().map(r -> r.getAttributeAsString(DuplicateSet.FGBIO_MOLECULAR_IDENTIFIER_TAG))
                .distinct().collect(Collectors.toList());
    }
}
