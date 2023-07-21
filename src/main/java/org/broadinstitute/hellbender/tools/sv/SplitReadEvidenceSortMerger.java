package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

import java.util.Comparator;

/**
 * Imposes additional ordering of same-locus SplitReadEvidence by sample and strand.
 * Imposes uniqueness criterion on <locus, sample, strand>.
 */
public class SplitReadEvidenceSortMerger extends AbstractEvidenceSortMerger<SplitReadEvidence> {

    public SplitReadEvidenceSortMerger( final SAMSequenceDictionary dictionary,
                                          final FeatureSink<SplitReadEvidence> outputSink ) {
        super(dictionary, outputSink,
                Comparator.comparing(SplitReadEvidence::getSample)
                    .thenComparing((f1, f2) -> Boolean.compare(f1.getStrand(), f2.getStrand())));
    }

    protected void complain( final SplitReadEvidence evidence ) {
        throw new UserException("Two instances of SplitReadEvidence for sample " +
                evidence.getSample() + " at " + evidence.getContig() + ":" +
                evidence.getStart() + (evidence.getStrand() ? " right" : " left"));
    }
}
