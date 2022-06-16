package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

import java.util.Comparator;

/**
 * Imposes additional ordering of same-locus LocusDepth records by sample.
 * Imposes uniqueness criterion on <locus, sample>.
 */
public class LocusDepthSortMerger extends AbstractEvidenceSortMerger<LocusDepth> {
    public LocusDepthSortMerger( final SAMSequenceDictionary dictionary,
                                  final FeatureSink<LocusDepth> outputSink ) {
        super(dictionary, outputSink, Comparator.comparing(LocusDepth::getSample));
    }

    protected void complain( final LocusDepth evidence ) {
        throw new UserException("Two instances of LocusDepth for sample " +
                evidence.getSample() + " at " + evidence.getContig() +
                ":" + evidence.getStart());
    }
}
