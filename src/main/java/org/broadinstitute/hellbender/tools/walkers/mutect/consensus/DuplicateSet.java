package org.broadinstitute.hellbender.tools.walkers.mutect.consensus;

import org.broadinstitute.hellbender.tools.walkers.mutect.UMI;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;


/** Maybe useful, maybe un**/
public class DuplicateSet {
    public static final String FGBIO_MOLECULAR_IDENTIFIER_TAG = "MI";
    public static final String FGBIO_MI_TAG_DELIMITER = "/";
    private int moleculeId = -1;
    private UMI umi;
    private String contig;
    private int fragmentStart = -1;
    private int fragmentEnd = -1;
    private List<GATKRead> reads;
    boolean smallInsert; // if true, the reads read into adaptors

    public DuplicateSet(){}

    public DuplicateSet(final GATKRead read){
        init(read);
    }

    public void init(GATKRead read){
        Utils.validate(moleculeId == -1 || moleculeId == getMoleculeID(read),
                String.format("Inconsisntent molecule IDs: Duplicate set id = %s, read molecule id = %s", moleculeId, getMoleculeID(read)));
        setMoleduleId(read);
        umi = new UMI(read);
        contig = read.getContig();
        if (read.isReverseStrand()){
            fragmentEnd = read.getEnd(); // TODO: does this include softclips?
        } else {
            fragmentStart = read.getStart();
        }
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

    private int getMoleculeID(final GATKRead read) {
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

    public SimpleInterval getDuplicateSetInterval(){
        Utils.validate(SimpleInterval.isValid(contig, fragmentStart, fragmentEnd), "Invalid duplicate set interval: " + new SimpleInterval(contig, fragmentStart, fragmentEnd));
        return new SimpleInterval(contig, fragmentStart, fragmentEnd);

    }

}
