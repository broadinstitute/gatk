package org.broadinstitute.hellbender.utils.fermi;

import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collection;

/** service provider interface for assemblers */
public interface Assembler {
    Assembly createAssembly( final Collection<GATKRead> reads );
}
