package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Locatable;

public interface Read extends Locatable {
    String getName();

    int getLength();

    byte[] getBases();

    byte[] getBaseQualities();

    Cigar getCigar();

    boolean isUnmapped();
}
