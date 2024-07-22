package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import org.broadinstitute.barclay.argparser.Argument;

import java.util.Comparator;
import java.util.Map;

public class CompareGtfInfo implements Comparator<Map.Entry<String, GtfInfo>>{

    @Argument()//TODO: figure this out
    SAMSequenceDictionary dictionary;

    CompareGtfInfo(SAMSequenceDictionary dictionary){
        this.dictionary = dictionary;
    }

//    @Override
//    public int compare(Map.Entry<String, GtfInfo> e1, Map.Entry<String, GtfInfo> e2) {
//        Interval e1Interval = e1.getValue().getInterval();
//        Interval e2Interval = e2.getValue().getInterval();
//        return Comparator.comparingInt((Interval interval) -> dictionary.getSequence(interval.getContig()).getSequenceIndex()).thenComparingInt(Interval::getStart).compare(e1Interval, e2Interval);
//    }


    @Override
    public int compare(Map.Entry<String, GtfInfo> e1, Map.Entry<String, GtfInfo> e2) {
        Interval e1Interval = e1.getValue().getInterval();
        Interval e2Interval = e2.getValue().getInterval();

        // Debugging: Print intervals and contig names
        //System.out.println("e1Interval: " + e1Interval);
        //System.out.println("e2Interval: " + e2Interval);

        String contig1 = e1Interval.getContig();
        String contig2 = e2Interval.getContig();

        //System.out.println("Contig1: " + contig1);
        //System.out.println("Contig2: " + contig2);

        // Get sequences
        SAMSequenceRecord sequenceRecord1 = dictionary.getSequence(contig1);
        SAMSequenceRecord sequenceRecord2 = dictionary.getSequence(contig2);

        //System.out.println("SequenceRecord1: " + sequenceRecord1);
        //System.out.println("SequenceRecord2: " + sequenceRecord2);

        // Get sequence indices
        Integer seqIndex1 = sequenceRecord1 != null ? sequenceRecord1.getSequenceIndex() : null;
        Integer seqIndex2 = sequenceRecord2 != null ? sequenceRecord2.getSequenceIndex() : null;

        //System.out.println("SeqIndex1: " + seqIndex1);
        //System.out.println("SeqIndex2: " + seqIndex2);

        // Use the sequence index in comparison
        return Comparator.comparingInt((Interval interval) -> {
                    SAMSequenceRecord record = dictionary.getSequence(interval.getContig());
                    if (record == null) {
                        System.out.println("No sequence found for contig: " + interval.getContig());
                        return Integer.MAX_VALUE; // Handle missing sequence case
                    }
                    return record.getSequenceIndex();
                }).thenComparingInt(Interval::getStart)
                .compare(e1Interval, e2Interval);
    }
}
