package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.Argument;

import java.util.Set;

/**
 * Keep only reads for a given sample.
 * Matching is done by case-sensitive exact match.
 */
public final class SampleReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = "sample_to_keep", shortName = "goodSM", doc="The name of the sample(s) to keep, filtering out all others", optional=false)
    public Set<String> samplesToKeep = null;

    @Override
    public boolean test(final SAMRecord read) {
        final SAMReadGroupRecord readGroup = read.getReadGroup();
        return readGroup != null && samplesToKeep.contains(readGroup.getSample());
    }
}
