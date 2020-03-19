package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.tools.walkers.mutect.consensus.DuplicateSet;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.stream.Collectors;

public abstract class DuplicateSetWalker extends ReadWalker {
    DuplicateSet currentDuplicateSet; // TODO: does it make sense to call clear() and recycle the same object?
    public static int INITIAL_MOLECULAR_ID = -1; // this won't match the ID of the first read

    @Override
    public final void onStartup(){
        // TODO: I had to make onStartUp of ReadWalker not final to do this. The right thing to do here is to
        // write a ReadWalkerBase parent class for both ReadWalker and DuplicateSet Walker, as Louis suggested.
        super.onStartup();
        currentDuplicateSet = new DuplicateSet();
        currentDuplicateSet.setMoleduleId(INITIAL_MOLECULAR_ID);
        // TODO: I guess the problem isn't really peeking but rather peeking gave me three reads instead of one...
    }

    /***
     * FGBio GroupByUMI returns reads sorted by molecular ID: For example, the input bam may look like
     * read1: ... MI:Z:0/A ...
     * read2: ... MI:Z:0/A ...
     * read3: ... MI:Z:0/B ...
     * read4: ... MI:Z:0/B ...
     * read5: ... MI:Z:1/A ...
     * read6: ... MI:Z:1/B ...
     * read7: ... MI:Z:1/B ...
     *
     * Thus it's sufficient to go through the reads in order and collect them in a list until
     * we encounter the next molecular ID, at which point we pass the list to the {@code apply}
     * method of the child class and clear the {@code currentDuplicateSet} variable.
     *
     */
    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (currentDuplicateSet.getMoleculeId() != DuplicateSet.getMoleculeID(read)) {
            if (rejectDuplicateSet(currentDuplicateSet)){
                currentDuplicateSet = new DuplicateSet(read);
                return;
            }

            apply(currentDuplicateSet,
                    new ReferenceContext(reference, currentDuplicateSet.getDuplicateSetInterval()), // Will create an empty ReferenceContext if reference or readInterval == null
                    new FeatureContext(features, currentDuplicateSet.getDuplicateSetInterval()));
            currentDuplicateSet = new DuplicateSet(read);
        } else {
            currentDuplicateSet.addRead(read);
        }
    }

    public abstract void apply(DuplicateSet duplicateSet, ReferenceContext referenceContext, FeatureContext featureContext );

    protected boolean rejectDuplicateSet(final DuplicateSet duplicateSet){
        if (!duplicateSet.hasValidInterval()){
            logger.info("Duplicate Set with Invalid Intervals");
            logger.info("Number of reads:" + currentDuplicateSet.getReads().size());
            if (currentDuplicateSet.getReads().size() > 0){
                logger.info("First read: " + currentDuplicateSet.getReads().get(0));
            }
            return true;
        } else {
            return false;
        }
    }

    @Override
    public void postProcess(){
        if (currentDuplicateSet.getReads().size() > 0){
            apply(currentDuplicateSet,
                    new ReferenceContext(reference, currentDuplicateSet.getDuplicateSetInterval()),
                    new FeatureContext(features, currentDuplicateSet.getDuplicateSetInterval()));
        }
    }
}
