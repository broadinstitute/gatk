package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.utils.read.GATKRead;

public interface ReadElementExtractor {
    String header();

    String extractElement(GATKRead read);
}

