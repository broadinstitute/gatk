package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

import java.util.Comparator;

/**
 * Imposes additional ordering of same-locus SiteDepth records by sample.
 * Imposes uniqueness criterion on <locus, sample>.
 */
public class SiteDepthSortMerger extends AbstractEvidenceSortMerger<SiteDepth> {
    public SiteDepthSortMerger( final SAMSequenceDictionary dictionary,
                                final FeatureSink<SiteDepth> outputSink ) {
        super(dictionary, outputSink, Comparator.comparing(SiteDepth::getSample));
    }

    protected void complain( final SiteDepth evidence ) {
        throw new UserException("Two instances of SiteDepth for sample " +
                evidence.getSample() + " at " + evidence.getContig() +
                ":" + evidence.getStart());
    }
}
