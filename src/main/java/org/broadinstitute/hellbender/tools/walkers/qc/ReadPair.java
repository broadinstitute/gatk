package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Data structure that contains the set of reads sharing the same queryname, including
 * the primary, secondary (i.e. multi-mapping), and supplementary (e.g. chimeric) alignments.
 *
 */
public class ReadPair {
    private GATKRead firstOfPair;
    private GATKRead secondOfPair;
    private final List<GATKRead> secondaryAlignments = new ArrayList<>();
    private final List<GATKRead> supplementaryAlignments = new ArrayList<>();
    private final String queryName;

    public ReadPair(final GATKRead read) {
        this.queryName = read.getName();
        add(read);
    }

    public int numberOfAlignments(){
        int num = 0;
        num = firstOfPair != null ? num + 1 : num;
        num = secondOfPair != null ? num + 1 : num;
        num += secondaryAlignments.size();
        num += supplementaryAlignments.size();
        return num;
    }

    public String getName() {
        return queryName;
    }

    public GATKRead getFirstOfPair(){
        return firstOfPair;
    }

    public GATKRead getSecondOfPair(){
        return secondOfPair;
    }

    public List<GATKRead> getSecondaryAlignments() {
        return secondaryAlignments;
    }

    public void add(final GATKRead read) {
        if (! this.queryName.equals(read.getName())){
            throw new UserException("Read names do not match: " + this.queryName + " vs " + read.getName());
        }

        if (isPrimaryAlignment(read) && read.isFirstOfPair()) {
            Utils.validate(this.firstOfPair != null,
                    "The primary firstOfPair is already set. Added read = " + read.getName());
            this.firstOfPair = read;
        } else if (isPrimaryAlignment(read) && read.isSecondOfPair()) {
            this.secondOfPair = read;
        } else if (read.isSecondaryAlignment()) {
            this.secondaryAlignments.add(read);
        } else if (read.isSupplementaryAlignment()) {
            this.supplementaryAlignments.add(read);
        } else {
            throw new UserException("Unknown read type: " + read.getContig());
        }
    }

    private boolean isPrimaryAlignment(final GATKRead read){
        return ! read.isSecondaryAlignment() && ! read.isSupplementaryAlignment();
    }

    public boolean isDuplicateMarked() {
        if (firstOfPair.isDuplicate()) {
            // Make sure the rest is duplicate-marked
            if (!secondOfPair.isDuplicate() || secondaryAlignments.stream().anyMatch(r -> !r.isDuplicate())) {
                throw new UserException("First of pair a duplicate but the rest is not" + secondOfPair.getName());
            }
        } else {
            // Make sure the rest is not duplicate-marked
            if (secondOfPair.isDuplicate() || secondaryAlignments.stream().anyMatch(r -> r.isDuplicate())) {
                throw new UserException("First of pair a not duplicate but the rest is " + secondOfPair.getName());
            }
        }

        return firstOfPair.isDuplicate();
    }

    public List<GATKRead> getReads(final boolean onlyPrimaryAlignments){
        List<GATKRead> reads = Arrays.asList(firstOfPair, secondOfPair);
        if (onlyPrimaryAlignments){
            return(reads);
        } else {
            reads.addAll(secondaryAlignments);
            reads.addAll(supplementaryAlignments);
            return(reads);
        }
    }
}
