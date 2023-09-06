package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

import java.util.Comparator;
import java.util.PriorityQueue;

/**
 * A FeatureSink that buffers and resolves (by merging, or by checking for redundancy) features
 * that occur on the same interval.
 */
public abstract class AbstractEvidenceSortMerger <F extends SVFeature> implements FeatureSink<F> {
    protected final SAMSequenceDictionary dictionary;
    protected final FeatureSink<F> outputSink;
    protected final PriorityQueue<F> sameLocusQueue;
    protected F currentLocus;

    public AbstractEvidenceSortMerger( final SAMSequenceDictionary dictionary,
                                       final FeatureSink<F> outputSink,
                                       final Comparator<? super F> comparator ) {
        this.dictionary = dictionary;
        this.outputSink = outputSink;
        this.sameLocusQueue = new PriorityQueue<>(comparator);
        this.currentLocus = null;
    }

    @Override
    public void write( final F feature ) {
        if ( currentLocus == null ) {
            currentLocus = feature;
            sameLocusQueue.add(feature);
        } else {
            final int cmp = IntervalUtils.compareLocatables(currentLocus, feature, dictionary);
            if ( cmp == 0 ) {
                sameLocusQueue.add(feature);
            } else if ( cmp < 0 ) {
                resolveSameLocusFeatures();
                currentLocus = feature;
                sameLocusQueue.add(feature);
            } else {
                throw new GATKException("features not presented in dictionary order");
            }
        }
    }

    @Override
    public void close() {
        resolveSameLocusFeatures();
        outputSink.close();
    }

    protected void resolveSameLocusFeatures() {
        if ( sameLocusQueue.isEmpty() ) {
            return;
        }
        final Comparator<? super F> comparator = sameLocusQueue.comparator();
        F lastEvidence = sameLocusQueue.poll();
        while ( !sameLocusQueue.isEmpty() ) {
            final F evidence = sameLocusQueue.poll();
            if ( comparator.compare(lastEvidence, evidence) == 0 ) {
                complain(evidence);
            }
            outputSink.write(lastEvidence);
            lastEvidence = evidence;
        }
        outputSink.write(lastEvidence);
    }

    protected abstract void complain( F feature );
}
