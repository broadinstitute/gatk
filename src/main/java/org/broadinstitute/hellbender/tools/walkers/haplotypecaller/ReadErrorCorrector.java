package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collection;
import java.util.List;

public interface ReadErrorCorrector {
    List<GATKRead> correctReads(final Collection<GATKRead> reads);
}
