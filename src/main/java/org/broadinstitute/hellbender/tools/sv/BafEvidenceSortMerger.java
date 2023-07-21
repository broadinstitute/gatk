package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

import java.util.Comparator;

/**
 * Imposes additional ordering of same-locus BafEvidence by sample.
 * Imposes uniqueness criterion on <locus, sample>.
 */
public class BafEvidenceSortMerger extends AbstractEvidenceSortMerger<BafEvidence> {
    public BafEvidenceSortMerger( final SAMSequenceDictionary dictionary,
                                  final FeatureSink<BafEvidence> outputSink ) {
        super(dictionary, outputSink, Comparator.comparing(BafEvidence::getSample));
    }

    protected void complain( final BafEvidence evidence ) {
        throw new UserException("Two instances of BafEvidence for sample " +
                                evidence.getSample() + " at " + evidence.getContig() +
                                ":" + evidence.getStart());
    }
}
