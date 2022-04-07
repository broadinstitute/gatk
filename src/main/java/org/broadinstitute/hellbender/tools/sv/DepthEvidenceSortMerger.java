package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

import static org.broadinstitute.hellbender.tools.sv.DepthEvidence.MISSING_DATA;

/**
 * Merges records for the same interval into a single record, when possible, throws if not possible.
 * It's assumed that all the records refer to the same samples in the same order.  (This can be
 * arranged by calling extractSamples on each record.)
 */
public class DepthEvidenceSortMerger implements FeatureSink<DepthEvidence> {
    private final SAMSequenceDictionary dictionary;
    private final FeatureSink<DepthEvidence> outputSink;
    private DepthEvidence mergedEvidence;

    public DepthEvidenceSortMerger( final SAMSequenceDictionary dictionary,
                                    final FeatureSink<DepthEvidence> outputSink ) {
        this.dictionary = dictionary;
        this.outputSink = outputSink;
        this.mergedEvidence = null;
    }

    @Override
    public void write( final DepthEvidence feature ) {
        if ( mergedEvidence == null ) {
            mergedEvidence = feature;
            return;
        }
        int cmp = IntervalUtils.compareLocatables(mergedEvidence, feature, dictionary);
        if ( cmp == 0 ) {
            merge(feature);
        } else if ( cmp < 0 ) {
            outputSink.write(mergedEvidence);
            mergedEvidence = feature;
        } else {
            throw new GATKException("features not presented in dictionary order");
        }
    }

    @Override
    public void close() {
        if ( mergedEvidence != null ) {
            outputSink.write(mergedEvidence);
            mergedEvidence = null;
        }
        outputSink.close();
    }

    private void merge( final DepthEvidence evidence ) {
        final int[] mergedCounts = mergedEvidence.getCounts();
        final int nCounts = mergedCounts.length;
        final int[] evidenceCounts = evidence.getCounts();
        if ( nCounts != evidenceCounts.length ) {
            throw new GATKException("All DepthEvidence ought to have the same sample list at this point.");
        }
        for ( int idx = 0; idx != nCounts; ++idx ) {
            final int count = evidenceCounts[idx];
            if ( count != MISSING_DATA ) {
                if ( mergedCounts[idx] == MISSING_DATA ) {
                    mergedCounts[idx] = count;
                } else {
                    throw new UserException("Multiple sources for count of sample#" + (idx+1) +
                            " at " + evidence.getContig() + ":" + evidence.getStart() + "-" +
                            evidence.getEnd());
                }
            }
        }
    }
}
