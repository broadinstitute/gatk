package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ReadPair {
    private GATKRead firstOfPair;
    private GATKRead secondOfPair;
    private List<GATKRead> secondaryAlignments = new ArrayList<>(10);
    private List<GATKRead> supplementaryAlignments = new ArrayList<>(10); // Finally understand the difference
    private String queryName = null;


    public ReadPair() { }

    public ReadPair(final GATKRead read) {
        this.queryName = read.getName();
        add(read);
    }

    public int size(){
        int size = 0;
        size = firstOfPair != null ? size + 1 : size;
        size = secondOfPair != null ? size + 1 : size;
        size += secondaryAlignments.size();
        return size += supplementaryAlignments.size(); // what should we do with the supplementary alignments?
    }

    public String getQueryName() {
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
        if (this.queryName == null){
            this.queryName = read.getName();
        }

        if (! this.queryName.equals(read.getName())){
            throw new UserException("Read names do not match: " + this.queryName + " vs " + read.getName());
        }

        if (isPrimaryAlignment(read) && read.isFirstOfPair()) {
            this.firstOfPair = read;
        } else if (isPrimaryAlignment(read) && read.isSecondOfPair()) {
            this.secondOfPair = read;
        } else if (read.isSecondaryAlignment()) {
            this.secondaryAlignments.add(read);
        } else if (read.isSupplementaryAlignment()) {
            this.supplementaryAlignments.add(read);
        } else {
            throw new UserException("Unknown read type");
        }
    }

    private boolean isPrimaryAlignment(final GATKRead read){
        return ! read.isSecondaryAlignment() && ! read.isSupplementaryAlignment();
    }

    public boolean isDuplicateMarked() {
        // Doing some investigation
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
