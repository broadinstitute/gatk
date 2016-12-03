package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;
import java.util.HashSet;
import java.util.Set;

/**
 * Keep only reads for a given sample.
 * Matching is done by case-sensitive exact match.
 */
public final class SampleReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = "sample",
            shortName = "sample", doc="The name of the sample(s) to keep, filtering out all others", optional=false)
    public Set<String> samplesToKeep = new HashSet<>();

    public SampleReadFilter( ) { }

    public SampleReadFilter( final SAMFileHeader header) {
        super.setHeader(header);
    }

    @Override
    public boolean test( final GATKRead read ) {
        final String sample = ReadUtils.getSampleName(read, samHeader);
        return sample != null && samplesToKeep.contains(sample);
    }
}
