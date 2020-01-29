package org.broadinstitute.hellbender.tools.walkers.mutect.consensus;

import org.broadinstitute.hellbender.tools.walkers.mutect.UMI;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;


/** Maybe useful, maybe un**/
public class DuplicateSet {
    private UMI umi;
    private int fragmentStart = -1;
    private int fragmentEnd = -1;
    private List<GATKRead> reads;
    boolean smallInsert; // if true, the reads read into adaptors

    public DuplicateSet(final GATKRead read){
        umi = new UMI(read);
        if (read.isReverseStrand()){
            fragmentEnd = read.getEnd(); // TODO: does this include softclips?
        } else {
            fragmentStart = read.getStart();
        }
    }

    public List<GATKRead> getReads(){
        return reads;
    }

    public boolean checkMembership(final GATKRead read){
        // check 1: do the UMI match
        final UMI readUMI = new UMI(read);
        if (!umi.equalsModuloOrder(readUMI)){
            return false;
        }

        // check 2: do they come from the same fragment?
        if (read.isReverseStrand()){
            if (fragmentEnd == -1){
                // START HERE, what kind of assumptions should we make on the input from groupByUMI here?
            }

        } else { // forward read
            return false; // placeholder
        }

        return true; // placeholder
    }

    public void addRead(final GATKRead read){
        return;
    }

}
