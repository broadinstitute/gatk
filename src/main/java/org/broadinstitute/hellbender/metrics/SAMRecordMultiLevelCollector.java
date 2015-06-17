package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.reference.ReferenceSequence;

/** Defines a MultilevelPerRecordCollector using the argument type of SAMRecord so that this doesn't have to be redefined for each subclass of MultilevelPerRecordCollector */
public abstract class SAMRecordMultiLevelCollector<BEAN extends MetricBase,
        HKEY extends Comparable<HKEY>> extends MultiLevelCollector<BEAN, HKEY, SAMRecord> {

    @Override
    protected SAMRecord makeArg(SAMRecord samRec, final ReferenceSequence refSeq) {
        return samRec;
    }
}
