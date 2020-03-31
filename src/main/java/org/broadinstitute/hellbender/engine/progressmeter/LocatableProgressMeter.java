package org.broadinstitute.hellbender.engine.progressmeter;

import htsjdk.samtools.util.Locatable;

/**
 * A ProgressMeter which accepts Locatables
 */
public class LocatableProgressMeter extends ProgressMeter<Locatable> {

    public LocatableProgressMeter(){
        super();
    }

    public LocatableProgressMeter(final double secondsBetweenUpdates ){
        super( secondsBetweenUpdates );
    }

    /**
     * @return genomic location of the most recent record formatted for output to the logger
     */
    @Override
    protected String formatRecord(Locatable locatable) {
        return locatable != null ? locatable.getContig() + ":" + locatable.getStart() : "unmapped";
    }
}
