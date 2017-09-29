package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;

public interface SmithWatermanAlignment {
    Cigar getCigar();
    int getAlignmentOffset();
}
