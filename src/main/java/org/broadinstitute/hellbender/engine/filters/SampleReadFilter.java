package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.HashSet;
import java.util.Set;

/**
 * Keep only reads for a given sample.
 * Matching is done by case-sensitive exact match.
 */
public final class SampleReadFilter implements ReadFilter, CommandLineFilter {
    private static final long serialVersionUID = 1L;

    private final static String sampleArgName = "sample";

    @Argument(fullName = sampleArgName, shortName = "goodSM", doc="The name of the sample(s) to keep, filtering out all others", optional=true)
    public Set<String> samplesToKeep = new HashSet<>();

    private final SAMFileHeader header;

    public SampleReadFilter( final SAMFileHeader header) {
        this.header = header;
    }

    @Override
    public boolean test( final GATKRead read ) {
        final String sample = ReadUtils.getSampleName(read, header);
        return sample != null && samplesToKeep.contains(sample);
    }

    @Override
    public String validate() {
        String message = null;
        if (samplesToKeep.size() <= 0) {
            message = "requires one or more values for \"" + sampleArgName + "\"";
        }

        return message;
    }

}
