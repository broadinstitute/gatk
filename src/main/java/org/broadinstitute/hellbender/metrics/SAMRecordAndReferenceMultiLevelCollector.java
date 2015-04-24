package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.reference.ReferenceSequence;
@SuppressWarnings("rawtypes")
public abstract class SAMRecordAndReferenceMultiLevelCollector<BEAN extends MetricBase,
        HKEY extends Comparable> extends MultiLevelCollector<BEAN, HKEY, SAMRecordAndReference> {

        @Override
        protected SAMRecordAndReference makeArg(SAMRecord samRec, final ReferenceSequence refSeq) {
            return new SAMRecordAndReference(samRec, refSeq);
        }
}


