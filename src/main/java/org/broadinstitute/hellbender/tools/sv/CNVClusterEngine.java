package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collection;

public class CNVClusterEngine extends LocatableClusterEngine<SVCallRecordWithEvidence> {
    public CNVClusterEngine(final SAMSequenceDictionary dictionary) {
        super(dictionary, CLUSTERING_TYPE.SINGLE_LINKAGE);
    }

    @Override
    protected boolean clusterTogether(SVCallRecordWithEvidence a, SVCallRecordWithEvidence b) {
        return false;
    }

    @Override
    protected SimpleInterval getClusteringInterval(SVCallRecordWithEvidence item, SimpleInterval currentClusterInterval) {
        return null;
    }

    @Override
    protected SVCallRecordWithEvidence deduplicateIdenticalItems(Collection<SVCallRecordWithEvidence> items) {
        return null;
    }

    @Override
    protected boolean itemsAreIdentical(SVCallRecordWithEvidence a, SVCallRecordWithEvidence b) {
        return false;
    }

    @Override
    protected SVCallRecordWithEvidence flattenCluster(Collection<SVCallRecordWithEvidence> cluster) {
        return null;
    }
}
