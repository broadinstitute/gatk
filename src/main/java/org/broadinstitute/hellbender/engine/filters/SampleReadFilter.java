package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;
import java.util.HashSet;
import java.util.Set;

/**
 * Keep only reads for a given sample.
 *
 * <p>Matching is done by case-sensitive exact match.</p>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only reads for a given sample")
public final class SampleReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = ReadFilterArgumentDefinitions.SAMPLE_NAME,
            shortName = ReadFilterArgumentDefinitions.SAMPLE_NAME, doc="The name of the sample(s) to keep, filtering out all others", optional=false)
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
