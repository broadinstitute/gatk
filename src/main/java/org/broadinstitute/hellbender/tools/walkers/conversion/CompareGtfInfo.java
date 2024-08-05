package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Comparator;
import java.util.Map;

public class CompareGtfInfo implements Comparator<Map.Entry<String, GtfInfo>> {

    @Argument(shortName = "SD", fullName = "SEQUENCE_DICTIONARY", optional =  true)//TODO: figure this out
    SAMSequenceDictionary dictionary;

    CompareGtfInfo(SAMSequenceDictionary dictionary) {
        this.dictionary = dictionary;
    }

    // compare two entries of a map where key = geneId or transcriptId and value = gtfInfo object
    @Override
    public int compare(Map.Entry<String, GtfInfo> e1, Map.Entry<String, GtfInfo> e2) {
        Interval e1Interval = e1.getValue().getInterval();
        Interval e2Interval = e2.getValue().getInterval();

        Utils.nonNull(dictionary.getSequence(e1Interval.getContig()), "could not get sequence for " + e1Interval.getContig());
        Utils.nonNull(dictionary.getSequence(e2Interval.getContig()), "could not get sequence for " + e2Interval.getContig());

        // compare by contig
        int contigComparison = Integer.compare(
                dictionary.getSequence(e1Interval.getContig()).getSequenceIndex(),
                dictionary.getSequence(e2Interval.getContig()).getSequenceIndex()
        );

        if (contigComparison != 0) {
            return contigComparison;
        }

        // compare by start if contigs are the same
        int startComparison = Integer.compare(e1Interval.getStart(), e2Interval.getStart());
        if (startComparison != 0) {
            return startComparison;
        }

        // tiebreaker: compare by key
        return e1.getKey().compareTo(e2.getKey());
    }
}
