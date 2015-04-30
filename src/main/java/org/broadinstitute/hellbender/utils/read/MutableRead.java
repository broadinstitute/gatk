package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Locatable;

public interface MutableRead extends Read {
    void setName( String name );

    void setBases( byte[] bases );

    void setBaseQualities( byte[] baseQualities );

    void setPosition( String contig, int start );

    void setPosition( Locatable locatable );

    void setCigar( Cigar cigar );

    void setCigar( String cigarString );

    void setIsUnmapped();
}