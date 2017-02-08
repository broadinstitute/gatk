package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;
import java.util.Comparator;

public class SimpleIntervalComparator implements Comparator<SimpleInterval>, Serializable {
    private static final long serialVersionUID = 1L;

    private final SAMSequenceDictionary referenceSequenceDictionary;

    public SimpleIntervalComparator(final SAMSequenceDictionary referenceSequenceDictionary) {
        this.referenceSequenceDictionary = referenceSequenceDictionary;
    }

    @Override
    public int compare(final SimpleInterval o1, final SimpleInterval o2) {
        final int contigComparison = new Integer(referenceSequenceDictionary.getSequenceIndex(o1.getContig())).compareTo(referenceSequenceDictionary.getSequenceIndex(o2.getContig()));
        if (contigComparison != 0) {
            return contigComparison;
        } else {
            return new Integer(o1.getStart()).compareTo(o2.getStart());
        }
    }
}
