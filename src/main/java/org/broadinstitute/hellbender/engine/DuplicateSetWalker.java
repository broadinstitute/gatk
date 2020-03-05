package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.tools.walkers.mutect.consensus.DuplicateSet;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public abstract class DuplicateSetWalker extends ReadWalker {
    DuplicateSet currentDuplicateSet; // TODO: does it make sense to call clear() and recycle the same object?
    public static int INITIAL_MOLECULAR_ID = -1; // this won't match the ID of the first read

    @Override
    public final void onStartup(){
        // TODO: I had to make onStartUp of ReadWalker not final to do this. The right thing to do here is to
        // write a ReadWalkerBase parent class for both ReadWalker and DuplicateSet Walker, as Louis suggested.
        super.onStartup();
        currentDuplicateSet = new DuplicateSet(peekFirstRead());
        currentDuplicateSet.setMoleduleId(INITIAL_MOLECULAR_ID);
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        // TODO: extract this logic and write a Duplicate Group Walker.
        if (!currentDuplicateSet.sameMolecule(read)) {
            // reference context is the reference context of the current read, not fragment, so it's note quite correct
            // to give it to letsDoIt
            // TODO: update the reference context to be the union of the span covered by the reads
            if (!currentDuplicateSet.hasValidInterval()){
                logger.info("Duplicate Set with Invalid Intervals");
                logger.info("Number of reads:" + currentDuplicateSet.getReads().size());
                if (currentDuplicateSet.getReads().size() > 0){
                    logger.info("First read: " + currentDuplicateSet.getReads().get(0));
                }
                logger.info("Intervals: " + currentDuplicateSet.getDuplicateSetInterval());
                // TODO: maybe collect statistics on discarded reads?
                // TODO: In the event of a really nasty read e.g. Cigar 25I121S...what should we do?
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

}
