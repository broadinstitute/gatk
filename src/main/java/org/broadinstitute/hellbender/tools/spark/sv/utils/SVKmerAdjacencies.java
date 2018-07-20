package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class SVKmerAdjacencies extends SVKmerLong {
    private int adjacentKmers;
    private static int[] reverseComplementAdjacencies;
    private static int[] bitCount = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };

    public SVKmerAdjacencies( final SVKmerLong kmer, final Base predecessor, final Base successor ) {
        super(kmer);
        adjacentKmers = (predecessor == null ? 0 : (0x10 << predecessor.ordinal())) |
                        (successor == null ? 0 : (1 << successor.ordinal()));
    }

    private SVKmerAdjacencies( final SVKmerLong kmer, final int adjacentKmers ) {
        super(kmer);
        this.adjacentKmers = adjacentKmers;
    }

    public int getAdjacencies() { return adjacentKmers; }

    @Override
    public SVKmerAdjacencies reverseComplement( final int kSize ) {
        return new SVKmerAdjacencies(super.reverseComplement(kSize), reverseComplementAdjacencies(adjacentKmers));
    }

    public void mergeAdjacencies( final SVKmerAdjacencies that ) {
        Utils.validateArg(equals(that), "merging unequal kmers");
        adjacentKmers |= that.adjacentKmers;
    }

    public int predecessorCount() { return bitCount[adjacentKmers >> 4]; }

    public int successorCount() { return bitCount[adjacentKmers & 0x0f]; }

    public boolean hasPredecessor( final Base base ) {
        return ((0x10 << base.ordinal()) & adjacentKmers) != 0;
    }

    public boolean hasSuccessor( final Base base ) {
        return ((1 << base.ordinal()) & adjacentKmers) != 0;
    }

    public SVKmerLong getSolePredecessor( final int kSize ) {
        switch (adjacentKmers >> 4) {
            case 1: return predecessor(Base.A, kSize);
            case 2: return predecessor(Base.C, kSize);
            case 4: return predecessor(Base.G, kSize);
            case 8: return predecessor(Base.T, kSize);
            default: return null;
        }
    }

    public SVKmerLong getSoleSuccessor( final int kSize ) {
        switch (adjacentKmers & 0x0f) {
            case 1: return successor(Base.A, kSize);
            case 2: return successor(Base.C, kSize);
            case 4: return successor(Base.G, kSize);
            case 8: return successor(Base.T, kSize);
            default: return null;
        }
    }

    public List<Integer> getPredecessorContigs( Map<SVKmerLong, Integer> contigEnds, final int kSize ) {
        final int nPredecessors = predecessorCount();
        if (nPredecessors == 0 ) return Collections.emptyList();
        List<Integer> predecessorList = new ArrayList<>(nPredecessors);
        final int predecessorBits = adjacentKmers >> 4;
        for ( int baseIdx = 0; baseIdx != 4; ++baseIdx ) {
            if ( (predecessorBits & (1 << baseIdx)) != 0 ) {
                final Integer contigId =
                        contigEnds.get(predecessor(baseValues[baseIdx], kSize).reverseComplement(kSize));
                if ( contigId == null ) {
                    throw new GATKException("can't find predecessor contig");
                }
                predecessorList.add(contigId);
            }
        }
        return predecessorList;
    }

    public List<Integer> getSuccessorContigs( Map<SVKmerLong, Integer> contigEnds, final int kSize ) {
        final int nSuccessors = successorCount();
        if ( nSuccessors == 0 ) return Collections.emptyList();
        List<Integer> successorList = new ArrayList<>(nSuccessors);
        final int successorBits = adjacentKmers & 0x0f;
        for ( int baseIdx = 0; baseIdx != 4; ++baseIdx ) {
            if ( (successorBits & (1 << baseIdx)) != 0 ) {
                final Integer contigId = contigEnds.get(successor(baseValues[baseIdx], kSize));
                if ( contigId == null ) {
                    throw new GATKException("can't find successor contig");
                }
                successorList.add(contigId);
            }
        }
        return successorList;
    }

    private static int reverseComplementAdjacencies( final int adjacentKmers ) {
        if ( reverseComplementAdjacencies == null ) {
            reverseComplementAdjacencies = new int[256];
            for ( int idx = 0; idx != 256; ++idx ) {
                int fwdVal = idx;
                int revVal = 0;
                for ( int bit = 0; bit != 8; ++bit ) {
                    revVal = (revVal << 1) | (fwdVal & 1);
                    fwdVal >>= 1;
                }
                reverseComplementAdjacencies[idx] = revVal;
            }
        }
        return reverseComplementAdjacencies[adjacentKmers];
    }
}
