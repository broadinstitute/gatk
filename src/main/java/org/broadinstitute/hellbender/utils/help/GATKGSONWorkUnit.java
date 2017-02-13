package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.barclay.help.GSONWorkUnit;

/**
 * Class representing a GSONWorkUnit for GATK work units.
 *
 * Adds "walkertype" to the base gson object created by Barclay.
 */
public class GATKGSONWorkUnit extends GSONWorkUnit {

    private String walkerType;

    public void setWalkerType(final String walkerType){
        this.walkerType = walkerType;
    }
}
