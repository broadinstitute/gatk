package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

public enum AssemblyRegionReadState {
    PRIMARY,        // This is the read's primary region
    NONPRIMARY,     // This region overlaps the read, but it is not primary
    EXTENDED,       // This region would overlap the read if it were extended
    UNMAPPED        // This read is not mapped
}