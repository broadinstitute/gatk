package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVFeature;
import org.broadinstitute.hellbender.utils.*;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Base class for querying SV evidence data.
 */

public abstract class SVEvidenceAggregator<T extends SVFeature> {

    private final FeatureDataSource<T> source;
    protected final SAMSequenceDictionary dictionary;

    /**
     * Constructor.
     *
     * @param source     feature source
     * @param dictionary sequence dictionary
     */
    public SVEvidenceAggregator(final FeatureDataSource<T> source,
                                final SAMSequenceDictionary dictionary) {
        Utils.nonNull(source);
        Utils.nonNull(dictionary);
        this.source = source;
        this.dictionary = dictionary;
    }

    /**
     * Returns the interval over which to query for evidence features for the given record. May return null.
     */
    abstract public SimpleInterval getEvidenceQueryInterval(final SVCallRecord record);

    /**
     * Returns true if the given evidence should be counted for the record. May assume that the evidence lies in
     * the interval returned by {@link #getEvidenceQueryInterval}.
     */
    abstract public boolean evidenceFilter(final SVCallRecord record, final T evidence);

    /**
     * Returns evidence features associated with the provided record.
     */
    public List<T> collectEvidence(final SVCallRecord call) {
        Utils.nonNull(call);
        final SimpleInterval callInterval = getEvidenceQueryInterval(call);
        if (callInterval == null) {
            return Collections.emptyList();
        }
        final List<T> rawEvidence = source.queryAndPrefetch(callInterval);

        // Filter overlapping evidence
        return rawEvidence.stream().filter(f -> evidenceFilter(call, f)).toList();
    }
}
