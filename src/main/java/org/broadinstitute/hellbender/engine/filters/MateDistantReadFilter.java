package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Keep only paired reads that are not near each other in a coordinate-sorted source of reads.
 * (I.e., they're both mapped, but mapped to different contigs, or are mapped too far apart on the same contig.)
 *
 * <p>See MappedReadFilter for criteria defining an mapped read.</p>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Keep only paired reads with mates mapped >= mate-too-distant-length (default 1KB) apart or on different contigs", extraDocs = ReadFilterLibrary.MappedReadFilter.class)
public class MateDistantReadFilter extends ReadFilter {
    public static final int DEFAULT_MATE_TOO_DISTANT_THRESHOLD = 1000;
    private static final long serialVersionUID = 1L;

    @Argument(fullName = ReadFilterArgumentDefinitions.MATE_TOO_DISTANT_LENGTH, doc="Minimum start location difference at which mapped mates are considered distant", optional=true)
    public int mateTooDistantLength = DEFAULT_MATE_TOO_DISTANT_THRESHOLD;

    public MateDistantReadFilter() {
    }

    public MateDistantReadFilter( final int minMappingQualityScore ) {
        this.mateTooDistantLength = minMappingQualityScore;
    }

    @Override public boolean test( final GATKRead read ) {
        return read.isPaired() && !read.isUnmapped() && !read.mateIsUnmapped() &&
                (Math.abs(read.getStart() - read.getMateStart()) >= mateTooDistantLength || !read.getContig().equals(read.getMateContig()));
    }
}
